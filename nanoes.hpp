#pragma once

#ifndef INCLUDED_NANOES_HPP
#define INCLUDED_NANOES_HPP

#include <cctype>
#include <string>
#include <vector>
#include <functional>
#include <map>
#include <tuple>
#include <algorithm>
#include <new>
#include <memory>
#include <numeric>
#include <assert.h>

#ifdef _MSC_VER
#include <malloc.h>
#else
#include <alloca.h>
#endif

//#define NANOES_GC_VERBOSE

#ifdef NANOES_GC_VERBOSE
#include <iostream>
#endif

namespace nanoes {
	class runtime;
	// we predeclare the precompiled class and the precompiler itself.
	class precompiled;
	class precompiler;

	// (see id template below) helper type identity templates used to do boxed type discovery
	template<class T>
	inline intptr_t id_all() {
		static int id;
		return intptr_t(&id);
	}
	// enforce specific id's for int/double
	template<>
	inline intptr_t id_all<int>() {
		return 1;
	}
	template<>
	inline intptr_t id_all<double>() {
		return 2;
	}

	// ref/const removing identity template.
	template<class T>
	inline intptr_t id() {
		return id_all<std::remove_const<std::remove_reference<T>::type>::type>();
	}

	// Values are stored as NAN tagged 64 bit ints
	// a value with all the high 32bits set is considered a reference instead of a double
	// where the lowest bit is a semi-space index and the rest is an offset into that space.
	// Use a union to indicate to the C compiler that we might access the values in 2 different ways.
	union INTERNALVALUE {
		uint64_t uval;
		double dval;
	};

	// nesvalue is a C++ root-reference that can be used from outside the runtime
	class nesvalue {
		// runtime pointer, used to identify what runtime the pointer belongs to.
		friend class runtime;
		friend class precompiled;
		runtime *rt;
		// prev-next links
		nesvalue *prev,*next;
		// the actual value
		INTERNALVALUE value;
		// regular linked list unlink/link functions
		void unlink() {
			if (prev) {
				prev->next = next;
				next->prev = prev;
				next = prev = nullptr;
			}
		}
		void link(nesvalue &other) {
			next = &other;
			prev = other.prev;
			prev->next = this;
			next->prev = this;
		}
		// private constructor only used by the inital runtime pointer.
		nesvalue(runtime *in_rt) : prev(this), next(this), rt(in_rt) {
			value.uval = 0;
		}
	public:
		// other refs may only be constructed from other root-refs
		nesvalue(const nesvalue& other) : rt(other.rt),value(other.value) {
			link(const_cast<nesvalue&>(other));
		}
		// the = operator connects root-refs to the right runtime.
		nesvalue& operator=(const nesvalue& other) {
			if (rt != other.rt) {
				unlink();
				link(const_cast<nesvalue&>(other));
				rt = other.rt;
			}
			value = other.value;
			return *this;
		}
		// TODO: operators for boxing/unboxing
		~nesvalue() {
			unlink();
		}
	};

	class runtime {
		// with them predeclared we can also friend them.
		friend class precompiled;
		friend class precompiler;

		// predeclare some classes used.
		class gcbase;
		class valuebase;
		template<class RT,class GCTX,class FTPL>
		struct parsectx;
		struct opbase;

		// the alignfter macro finds the first aligned pointer on or after the given one.
		static inline void* alignafter(size_t align, void *ip) {
			// find the next possible free position
			uintptr_t p = (uintptr_t)ip;
			// align our pointer properly.
			p = ((p + align - 1) / align)*align;
			return (void*)p;
		}

		// represent each semi-space with the space structure.
		struct space {
			size_t sz = 0;
			char * ptr=nullptr; // this is the real ptr, hptr is an ptr that might be offset a bit.
			gcbase *head = nullptr;  // first object in heap (might not be same as ptr due to alignment)
			gcbase *tail = nullptr;  // last active object in heap, used to figure out the next free location.
			bool owns(void *p) {
				uintptr_t start = (uintptr_t)ptr;
				uintptr_t end = start+sz;
				uintptr_t test = (uintptr_t)p;
				return start <= test && test < end;
			}
			// this function 
			void killall() {
				if (head) {
					bool last;
					do {
						gcbase *toDel = head; // object to delete
						head = head->next;       // next object
						toDel->~gcbase();     // delete it
						last = head == tail;
					} while (last);              // continue until we've done with the tail.
					head = tail = nullptr;
				}
			}
			inline void* fitalignafter(size_t dest_size, size_t align, void *ip) {
				uintptr_t p = (uintptr_t)alignafter(align, ip);
				// ensure that we won't overflow by using that address.
				if (p + dest_size > (uintptr_t)ptr + sz) {
					return nullptr;
				}
				return (void*)p;
			}
			// due to alignment and unknown objects we locate the next free position
			inline void* findfree(size_t dest_size,size_t align) {
				return fitalignafter(dest_size,align,tail ? ((char*)tail) + tail->size() : ptr);
			}
			template<class T>
			inline T* findfree() {
				return (T*)findfree(sizeof(T), alignof(T));
			}
			~space() {
				killall();
				if (ptr)
					delete[] ptr;
			}
		};
		space heap[2];                         // heap info (used for heap-ops)
		char *hptr[2] = { nullptr,nullptr };   // offset heap ptrs (these are used mostly by the unboxing/boxing operations)
		int aHeap = 0;                         // indicates the active heap


		// lexer/parser info and a complete symbol to key mapping table.
		std::vector<std::pair<int, std::string>> toks;
		int nextTok = 256;
		int atok(const char* v) {
			auto where = std::find_if(toks.begin(), toks.end(), [&](std::pair<int, std::string>& k) { return v == k.second; });
			if (where == toks.end()) {
				int out = nextTok++;
				toks.emplace_back(out, v);
				toksort();
				return out;
			} else {
				return where->first;
			}
		}

		// 1 : seq
		const int SEMICOLON = atok(";"), COMMA = atok(","), LBRA = atok("{"), RBRA = atok("}");
		// 2 : yield
		// 3 : assign
		// 4 : conditional
		// 5 : logor
		// 6 : logand
		// 7 : bitor
		// 8 : bitxor
		// 9 : bitand
		// 10 : equality
		const int EQ = atok("=="), NEQ = atok("!=");
		// 11 : relative comparison
		const int LT = atok("<"), LEQ = atok("<="), GEQ = atok(">="), GT = atok(">");
		// 13 : additive
		const int ADD = atok("+"), T_SUB = atok("-");
		// 14 : multiplicative
		const int MUL = atok("*"), T_DIV = atok("/"), T_MOD = atok("%");
		// 15 : exponentiation
		// 16 : prefix ops..
		// 17 : postdec
		// 18 : new (no-arg)
		// 19 : fncall/access
		const int LPAR = atok("("), RPAR = atok(")"), DOT=atok(".");
		
		const int T_IF = atok("if"), T_ELSE = atok("else"), T_WHILE=atok("while"), T_FUNCTION = atok("function"), T_RETURN = atok("return"), T_THIS = atok("this");
		const int ID = atok("   _ID_"), DNUM = atok("   _dnum_"), STR = atok("   __str__");
		int USER; // USER indicates the first user created token index!

		// the toks list MUST be sorted for the lexer to function
		void toksort() {
			std::sort(toks.begin(), toks.end(), [](auto a, auto b) {  return a.second < b.second; });
		}



		class gcbase;
		typedef runtime toucher;

		std::map<int, INTERNALVALUE> global;

		class gcbase {
			friend runtime;
			gcbase *next = nullptr;
			gcbase *moved = nullptr;
		public:
			virtual ~gcbase() {}
			virtual size_t size() = 0;
			virtual size_t align() = 0;
			virtual void moveto(void *p) = 0;
			virtual void touch(toucher& perfield) = 0;
		};

		class valuebase : public gcbase {
		protected:
			virtual intptr_t ty() = 0;
			valuebase() {}
		public:
			virtual void * get() = 0;

			template<class T> inline T* get() {
				if (ty() == id<T>())
					return (T*)get();
				else
					return nullptr;
			}
			virtual INTERNALVALUE invoke(runtime &rt,int argc, INTERNALVALUE * args) = 0;
			virtual std::string to_string() = 0;
			virtual double to_num() = 0;
			virtual bool truthy() = 0;

			// basic operations: prop-ref?
			// vtable? (we could need an vtable upgrade thing but would buy us a lot of speed?)
			// - VT value read: { let off = obj->vt[PIDX]; v = off?*(obj+off):obj->getter(PN) }
			// - VT value write { let off = obj->vt[PIDX]; off?(*(obj+off)=v):obj->setter(PN,v) }
			// Faster C++ integration:
			// - RPOCO2(-like?) auto-VT-gen
			//  - plus side of being able to keep WebGL refs as shared_ptr's that auto-unwraps for the native calls.
			// - manual possible?
		};

		inline valuebase* to_ptr(INTERNALVALUE v) {
			if (v.uval < (uint64_t)0xffffffff00000000LL) {
				return 0;
			} else {
				uint32_t i = (uint32_t)v.uval;
				return (valuebase*)(hptr[i & 1] + i);
			}
		}

		template<class T>
		inline T unbox(INTERNALVALUE v) {
			if (valuebase *p = to_ptr(v)) {
				if (p->ty() != id<T>()) {
					throw std::runtime_error(std::string("Invalid type in value"));
				}
				return *(T*)p->get();
			} else {
				throw std::runtime_error(std::string("Tried getting a ")+std::string(typeid(T).name())+" but got a number");
			}
		}
		template<>
		inline void* unbox<void*>(INTERNALVALUE v) {
			if (valuebase *p = to_ptr(v)) {
				return p;
			} else {
				return nullptr;
			}
		}
		template<>
		inline double unbox<double>(INTERNALVALUE v) {
			if (v.uval < (uint64_t)0xffffffff00000000LL) {
				return v.dval;
			} else {
				return to_ptr(v)->to_num();
			}
		}
		template<>
		inline std::string unbox<std::string>(INTERNALVALUE v) {
			if (valuebase *p = to_ptr(v)) {
				return p->to_string(); // TODO: do we want to move to GC strings?
			} else {
				return std::to_string(unbox<double>(v));
			}
		}
		template<>
		inline int unbox<int>(INTERNALVALUE v) {
			return (int)unbox<double>(v);
		}

		inline INTERNALVALUE invoke(INTERNALVALUE fn, int arg, INTERNALVALUE* args) {
			if (valuebase *p = to_ptr(fn)) {
				return p->invoke(*this, arg, args);
			} else {
				throw std::runtime_error(std::string("number ") + std::to_string(fn.dval) + " is not a function");
			}
		}

		inline bool truthy(INTERNALVALUE v) {
			if (v.uval == V_TRUE.uval) return true;
			if (v.uval == V_FALSE.uval) return false;
			if (v.uval == V_NULL.uval) return false;
			if (valuebase *p = to_ptr(v)) {
				return p->truthy();
			} else {
				return unbox<double>(v)!=0;
			}
		}

		// 1st pass, create scopes and note scopes. Also register usages.
		// mid-pass, find variable targets and mark scopes non-local.
		// 2nd pass, just use scope-info to 
		struct scopeinfo {
			bool top = false;
			std::shared_ptr<scopeinfo> parent = nullptr;
			std::vector<std::pair<int, std::string>> names;
			size_t maxnames=0;
			size_t maxstack = 0;
			bool local = true; // does variables in this scope escape? in that case we always need to allocate these scopes on the heap!
		};



		struct scope : public gcbase {
			virtual size_t size() {
				return (size_t)(((INTERNALVALUE*)alignafter(align(), (void*)sizeof(*this))) + info->maxnames);
			}
			size_t runtime::scope::align()
			{
				return std::max(alignof(decltype(*this)), alignof(INTERNALVALUE));
			}
			virtual void moveto(void* dest) {
				scope * out = (scope*)dest;
				INTERNALVALUE* slots = (INTERNALVALUE*)alignafter(align(), out);
				memcpy(slots, this->nslots, info->maxnames);
				::new (out) scope(slots, info);
			}
			virtual void touch(toucher& to) {
				to.touch(parent);
				int max = info->maxnames;
				for (int i = 0;i < max;i++)
					to.touch(nslots[i]);
			}

			scope * parent = nullptr;

			scopeinfo *info;

			INTERNALVALUE* nslots;

			scope(INTERNALVALUE* islots, scopeinfo *in_info) : nslots(islots), info(in_info) {
				if (islots)
					for (size_t i = 0;i < in_info->maxnames;i++)
						islots[i].uval = 0; // nunpat flips to zero
			}
		};
	
		struct frame;
		frame *activeframe= nullptr;
		struct frame {
			runtime *rt;
			frame *prev;
			scope** mscope;
			int ssize;
			INTERNALVALUE* stack;
			frame(runtime *in_rt, scope** in_scope,int issize,INTERNALVALUE* istack) : mscope(in_scope), rt(in_rt), prev(rt->activeframe),ssize(issize),stack(istack) {
				rt->activeframe = this;
				for (int i = 0;i < issize;i++)
					stack[i].uval = 0;
			}
			~frame() {
				assert(rt->activeframe == this);
				rt->activeframe = prev;
			}
		};


		INTERNALVALUE V_TRUE;  // TODO: set and collect
		INTERNALVALUE V_FALSE; // TODO: set and collect
		INTERNALVALUE V_NULL;  // TODO: set and collect
		nesvalue uroot{ this };

		template<class T>
		void touch(T* &v,bool allownew=false) {
			gcbase *old = v;
			if (!old)
				return; // no need updating null ptrs..
			// this object indicates that it's already old and "moved"
			if (old->moved) {
				v = (T*)old->moved;
				return;
			}
#ifdef NANOES_GC_VERBOSE
			std::cout << "Moving a " << typeid(T).name() << " from " << (void*)v << "  ";
#endif
			if (!heap[aHeap ^ 1].owns(old)) {
				if (allownew)
					return;
				throw std::exception("target not in old-gen");
			}
			if (heap[aHeap ^ 1].owns(&v))
				throw std::exception("old-loc is in old-gen");
			space *tospace = heap + aHeap;
			// find a new placement!
			auto targetpos=tospace->findfree(old->size(), old->align());
			if (!targetpos)
				throw std::exception("OUT OF MEMORY DURING GC!!!");
			// ask the object to move itself.
			old->moveto(targetpos);
			// update the passed pointer and the forwarding address.
			old->moved = v = (T*)targetpos;
			// link in the new object to the target-heap.
			if (tospace->tail) {
				tospace->tail->next = v;
			} else {
				tospace->head = v;
			}
			tospace->tail = v;
#ifdef NANOES_GC_VERBOSE
			std::cout << "was moved to " << (void*)v << std::endl;
#endif
		}
		void touch(INTERNALVALUE& vr,bool allownew=false) {
			// only need to work on pointers, so change to a pointer and then rebase it..
			if (auto *p = to_ptr(vr)) {
				touch(p,allownew);
				vr=ptr_to_val(p);
			}
		}
		void collect() {
#ifdef NANOES_GC_VERBOSE
			std::cout << "***** Begun collection" << std::endl;
#endif
			space *fromspace = heap + aHeap;
			// switch heaps and begun processing!!
			aHeap ^= 1;
			// mark some builtins
			touch(V_NULL);
			touch(V_FALSE);
			touch(V_TRUE);
			// mark the global objects
			for (auto& glob : global) {
				touch(glob.second);
			}
			// mark the user root set.
			nesvalue *croot = &uroot;
			do {
				touch(croot->value);
				croot = croot->next;
			} while (croot != &uroot);
			// mark the active stack frames and associated scopes.
			frame *t = activeframe;
			while (t) {
				for (int i = 0;i < t->ssize;i++) {
					touch(t->stack[i] );
				}
				scope* s = *t->mscope;
				if (s->info->local) { // local scope, only touch the currently active stack entries.
					s->touch(*this);  // delegate touching to the scope-obj itself
					assert(!s->parent || !s->parent->info->local); // ensure that a local cannot point to a parent-local
				} else {
					touch(*t->mscope);
				}
				t = t->prev;
			}

			space *tospace = heap + aHeap;
			gcbase * cur = tospace->head;
			assert(cur);
			while (cur) {
				cur->touch(*this);
				cur = cur->next;
			}

			// destroy everything left in the from-space
			fromspace->killall();

#ifdef NANOES_GC_VERBOSE
			int occ = ((uintptr_t)tospace->tail) - ((uintptr_t)tospace->ptr);
			std::cout << "***** Finished collection:";
			std::cout << "Occupacy:" << occ << " Tot:" << tospace->sz << "(" << (100.f*occ / (float)tospace->sz) <<")"<< std::endl;
#endif
		}

		template<class T>
		class value : public valuebase {
			T data;

			template<class T>
			void touch(T&,toucher& to) {
#ifdef NANOES_GC_VERBOSE
				std::cout<<"Ignoring to touch "<<typeid(T).name()<<std::endl;
#endif
			}
			
			template<class R, class ... ARGS, size_t... Is>
			INTERNALVALUE invoke_impl(runtime &rt, std::function<R(ARGS...)> & data, int argc, INTERNALVALUE * args, std::index_sequence<Is...>) {
				return rt.box(  data(   rt.unbox<ARGS>( args[Is] )...  ));
			}
			template<class R, class ... ARGS>
			INTERNALVALUE invoke(runtime &rt, std::function<R(ARGS...)> & data, int argc, INTERNALVALUE * args) {
				if (argc != sizeof...(ARGS)) {
					// TODO: location
					throw std::exception("Wrong script function arity");
				}
				return invoke_impl(rt,data, argc, args, std::make_index_sequence<sizeof...(ARGS)>{});
			}
			template<class T>
			INTERNALVALUE invoke(runtime &rt, T& ,int argc, INTERNALVALUE * args) {
				throw std::runtime_error(std::string("Trying to invoke on a ")+typeid(T).name());
			}

			double to_num(std::string & s) {
				throw std::runtime_error("Not a number : [" + s + "]");
			}
			template<class T>
			double to_num(T & v) {
				throw std::runtime_error(typeid(v).name() + std::string(" is not a number"));
			}

			bool truthy(std::string& r) {
				return !r.empty();
			}
			bool truthy(bool v) {
				return v;
			}
			template<class T>
			bool truthy(T & v) {
				return true;
			}

			std::string to_string(bool v) {
				return v ? "true" : "false";
			}
			std::string to_string(std::string & s) {
				return s;
			}
			template<class T>
			std::string to_string(T& arg) {
				return typeid(arg).name();
			}
		public:
			value() : data() {}
			value(const T & indata) : data(indata) {}
			value(T && indata) : data(indata) {}
			virtual ~value() {}
			virtual intptr_t ty() { return id<T>(); }
			virtual size_t size() { return sizeof(*this); }
			virtual size_t align() { return alignof(decltype(*this)); }
			virtual void touch(toucher& perfield) {
				touch(data, perfield);
			}
			virtual void *get() { return &data; }
			virtual void moveto(void* dest) {
#ifdef NANOES_GC_VERBOSE 
				std::cout << "[["<< typeid(*this).name() <<"]]";
#endif
				::new(dest) value(std::move(*this));
			}
			virtual INTERNALVALUE invoke(runtime &rt,int argc, INTERNALVALUE * args) {
				return invoke(rt,data, argc, args);
			}
			virtual double to_num() {
				return to_num(data);
			}
			virtual std::string to_string() {
				return to_string(data);
			}
			virtual bool truthy() {
				return truthy(data);
			}
		};

		struct prop_slot {
			INTERNALVALUE prop;
			INTERNALVALUE value;
		};
		template<int I>
		class dynobj : public valuebase {
			std::unique_ptr<std::vector<prop_slot>> extras;
		};


		INTERNALVALUE ptr_to_val(void *p) {
			uintptr_t ip = (uintptr_t)p;
			auto dist = ip - ((uintptr_t)hptr[aHeap]);
			if (dist >= heap[aHeap].sz) {
				abort(); // TODO: fix so that the ptr_to_val macro can detect ptrs from the other heap (although that should never be done/used!)
			}
			INTERNALVALUE out;
			out.uval = 0xffffffff00000000LL | ( dist|aHeap );
			return out;
		}


		template<class T>
		inline typename std::enable_if< std::is_base_of<valuebase, T>::value, INTERNALVALUE>::type box(T && iv) {
			T *p = heap[aHeap].findfree<T>();
			if (!p) {
				collect();
				p = heap[aHeap].findfree<T>();
				if (!p) {
					throw std::runtime_error(std::string("Out of GC heap while trying to allocate an ")+std::string(typeid(T).name()));
				}
			}
			// now construct the value in-place..
			new(p) T(std::move(iv));
			// and now update the head/tail for the heap.
			if (heap[aHeap].tail) {
				heap[aHeap].tail->next = p;
			} else {
				heap[aHeap].head = p;
			}
			heap[aHeap].tail = p;

			return ptr_to_val(p);
		}
		template<class T>
		inline typename std::enable_if< !std::is_base_of<valuebase, T>::value, INTERNALVALUE>::type box(T && iv) {
			return box(value<T>(std::move(iv)));
		}

		inline INTERNALVALUE box(double v) {
			INTERNALVALUE vb;
			vb.dval = v;
			return vb;
			//if (vb.uval < 0xffffffff00000000LL) {
			//	// if the pattern ok then return it
			//	return vb;
			//} else {
			//	// we had a NaN pattern that collided with our nunmask, return another nan
			//	vb.uval &= 0xfffff8ff00000000LL;
			//	return vb;
			//}
		}
		INTERNALVALUE box(int v) {
			return box((double)v);
		}

		enum class OPC : uint8_t {
			URETURN=0,
			RETURN=1,
			LOADINT=2,      // 2 value : next 8 is target-slot
			LOADDOUBLE=3,   // 3 values, 8 bits of target slot, 2 following values that are low/high bits
			LOADLIT=4,      // 1 value : 8bits op, 8 target, ????
			LOADGLOBAL=5,     // 2 values :8bits op, 8 target, 16bit target idx
			LOADNULL=6,     // 1 value  : 8bits op, 8 target
			LOADTRUE=7,
			LOADFALSE=8,
			LOADSSYM=9,       // 2 values :8bits op, 8bits target, 10 bits sym-idx and 6 bit up-count
			INVOKE=0xa,       // 1 value : 8bit ops, 8bits base, 16 bits count
			// 0xb-0xf
			ADD, SUB, MUL, DIV, MOD, // pops 2 and computes something. (8bit op, 8bit target)
			// 0x10-0x15
			LT, LEQ, GEQ, GT, EQ, NEQ,
			// 0x16
			FGOTO,        // 1 value : 8bit op, 8bit slot, 16bit jump offset
			// 0x17
			GOTO,         // 1 value : 8bit op, 8bit ????, 16bit jump offset
			LOADSPROP,    // "3" values: 8bit op, 8bit slot, 16bit ??  : 32bit prop id : per-load-ptr
			STORESPROP,   // "3" values: 8bit op, 8bit slot, 16bit ??  : 32bit prop id : per-store-ptr
		};

		class funtpl;
		struct fungenctx {
			runtime *rt;
			std::vector<int32_t> code;
			std::vector<nesvalue> literals;
			std::vector<int> globals;
			std::vector<std::shared_ptr<funtpl>> functions;
			size_t sp = 0;


			fungenctx(runtime *irt) : rt(irt) {}

			void reset() {
				sp = 0;
				literals.clear();
				globals.clear();
				code.clear();
				functions.clear();
			}

			void addfn(std::shared_ptr<funtpl> && fn) {
				functions.push_back(std::move(fn));
			}
			int addlit(const std::string& v) {
				int lid = literals.size();
				literals.emplace_back(rt->uroot);
				literals[lid].value = rt->box<std::string>(std::string(v));
				return lid;
			}
			void add(OPC a, int32_t b = 0, int32_t c = 0, bool reladdr = false, double *dvp = nullptr) {
				if (reladdr) {
					//printf("Goto label %d from %d\n",c,gctx->code.size());
					c -= code.size() + 1;
				} else {
					//printf("Op %d at %d\n", a, gctx->code.size());
				}
				int32_t v = (((int)a) | (b << 8) | (c << 16));
				code.push_back(v);
				if (dvp) {
					INTERNALVALUE v = rt->box(*dvp);
					code.push_back(int32_t(v.uval));
					code.push_back(int32_t(v.uval >> 32));
				}
			}
			int label() {
				//printf("Label at %d\n", gctx->code.size());
				return code.size();
			}
		};

		class funtpl {
			friend runtime;
			runtime *rt;
			std::pair<int, std::string> id;
			std::shared_ptr<scopeinfo> info;
			int argcount;
			std::vector<int32_t> code;
			std::vector<nesvalue> literals;
			std::vector<INTERNALVALUE*> globals;
			std::vector<std::shared_ptr<funtpl>> functions;
			std::string src;
		public:
			funtpl(runtime *in_rt, std::pair<int, std::string> in_id, int in_argcount, std::shared_ptr<scopeinfo>& in_si, fungenctx & in_gctx,std::string&& insrc)
				: rt(in_rt), id(in_id), argcount(in_argcount), info(in_si),code(std::move(in_gctx.code)),literals(in_gctx.literals),functions(std::move(in_gctx.functions)),src(insrc)
			{
				for (auto id : in_gctx.globals) {
					globals.push_back(&rt->global[id]);
				}
			}
		};
		struct funinst : valuebase {
			std::shared_ptr<funtpl> code;
			funtpl* qp;
			scope * parent;
			funinst(scope *in_parent, const std::shared_ptr<funtpl>& in_code) : parent(in_parent), code(in_code), qp(in_code.get()) {}
			virtual ~funinst() {}
			virtual size_t size() { return sizeof(*this); }
			virtual size_t align() { return alignof(decltype(*this)); }
			virtual void moveto(void *p) { ::new(p) funinst(std::move(*this)); }
			virtual void touch(toucher& to) {
				to.touch(parent);
			}
			virtual intptr_t ty() { return id<funinst>(); }
			virtual void * get() { return nullptr; } // no embedded objects here

			virtual std::string to_string() { return "<FUNCTION...>"; }
			virtual double to_num() { return 0; }
			virtual bool truthy() {return true;}
			static inline INTERNALVALUE invoke_impl(runtime *rt, funinst *fi, int argc, INTERNALVALUE *args);
			virtual INTERNALVALUE invoke(runtime &rt, int argc, INTERNALVALUE *args) {
				return invoke_impl(&rt, this, argc, args);
			}

		};

		template<class RT,class GCTX,class FTPL>
		struct parsectx {
			RT * rt;
			const char *state;
			std::string tok;

			parsectx(RT *inrt):rt(inrt) {}

			std::vector<std::pair<scopeinfo*, int>> accesses;
			std::vector<std::pair<const char*, std::shared_ptr<scopeinfo>>> scopes;
			std::shared_ptr<scopeinfo> topscope = nullptr;
			std::shared_ptr<scopeinfo> curscope = nullptr;
			bool prep = true;

			std::map<void*, int> cgmem;
			GCTX * gctx;

			void add(OPC a, int32_t b = 0, int32_t c = 0, bool reladdr = false, double *dvp = nullptr) {
				gctx->add(a, b, c, reladdr, dvp);
			}

			struct spres {
				parsectx *ctx;
				size_t off;
				size_t sz;
				spres(parsectx &ictx, int x) : ctx(&ictx), off(ctx->gctx->sp), sz(x) {
					ctx->gctx->sp += sz;
					scopeinfo * ascope = ctx->curscope ? ctx->curscope.get() : ctx->topscope.get();
					ascope->maxstack = std::max(ctx->gctx->sp, ascope->maxstack);

				}
				~spres() {
					ctx->gctx->sp -= sz;
					assert(ctx->gctx->sp == off);
				}
			};


			std::shared_ptr<scopeinfo> enter_scope() {
				auto loc = std::find_if(scopes.begin(), scopes.end(), [this](auto& v) { return v.first == state; });
				if (loc != scopes.end()) {
					return curscope = loc->second;
				} else {
					auto neu = std::make_shared<scopeinfo>();
					neu->parent = curscope;
					scopes.emplace_back(state, neu);
					return curscope = neu;
				}
			}
			void leave_scope() {
				curscope = curscope->parent;
			}
			std::pair<int, int> access(int tok, scopeinfo * from = nullptr, int count = 0) {
				if (!from)
					from = curscope.get();
				if (prep) {
					accesses.emplace_back(from, tok);
				} else if (from) {
					auto loc = std::find_if(from->names.begin(), from->names.end(), [&](auto& v) { return v.first == tok; });
					if (loc != from->names.end()) {
						if (count)
							from->local = false; // not a local access, taint it!
						return std::make_pair(count, int(loc - from->names.begin()));
					}
					if (scopeinfo * parent = from->parent.get()) {
						auto up = access(tok, parent, count + 1);
						if (up.first != -1 && count) {
							from->local = false; // variable found in sub, taint this one also if it was non-local.
						}
						return up;
					}
				}
				return std::make_pair(-1, -1);
			}

			// this function indicates if the character is a valid symbol character
			int sym(int x, bool first = true) { return x == '_' || (first ? isalpha(x) : isalnum(x)); }

			// eat one token from the stream into the context
			int lexeat(parsectx & ctx, double *ddest = nullptr) {
				const char *& state = ctx.state;
				// skip over spaces (TODO: update line numbers?)
				while (std::isspace(*state)) state++;
				// search our sorted token+ID list by narrowing
				int rs = 0, re = rt->toks.size() - 1;
				for (size_t sz = 0;rs <= re;) {
					// pre-token is smaller, go to next
					if (rt->toks[rs].second.size() != sz && (rt->toks[rs].second.size()<sz || rt->toks[rs].second.c_str()[sz] < state[sz])) {
						rs++;
						continue;
					}
					// post-token is larger, go to next
					if (rt->toks[re].second.size() != sz && (rt->toks[re].second.size() < sz || rt->toks[re].second.c_str()[sz] > state[sz])) {
						re--;
						continue;
					}
					// pre- and post- token are equal and we've found it. (also make sure we don't have an ID token and it continues!)
					if (rs == re && sz == rt->toks[rs].second.size() && !(sym(state[0]) && sym(state[sz], false))) {
						ctx.tok.assign(state, sz);
						state += sz;
						return rt->toks[rs].first;
					}
					// ok, we haven't found a token or needed to narrow our range to so add another character
					sz++;
				}
				// ok not an known operator or identifier, so either we have a new identifier or some kind of literal like a string or number
				ctx.tok.clear();
				if (*state == '\"') {
					// parse a strings
					state++;
					char c;
					while ((c = *state) != '\"') {
						if (!c)
							return -1;
						if (c == '\\') {
							state++;
							if (*state == 'n') c = '\n';
							else if (*state == 'r') c = '\r';
							else if (*state == 'b') c = '\b';
							else if (*state == 't') c = '\t';
							// TODO: more charcter escape codes. ?
							else c = *state;
						}
						ctx.tok.push_back(c);
						state++;
					}
					state++;
					return rt->STR;
				} else if (isdigit(*state) || (*state == '-'&&isdigit(state[1]))) {
					const char * start = state;
					double num = 0;
					double sign = *state == '-' ? -1 : 1;
					if (sign < 0)
						state++;
					while (isdigit(*state)) {
						num = (num * 10) + (*state - '0');
						state++;
					}
					ctx.tok.assign(start, state - start);
					if (ddest) {
						*ddest = num*sign;
					}
					return rt->DNUM;
				} else if (sym(*state)) {
					while (sym(*state, false))
						ctx.tok.push_back(*state++);
					int rv = rt->atok(ctx.tok.c_str());
					return rv;
				} else if (*state == 0) {
					ctx.tok = "<END OF FILE>";
					return -1;
				} else {
					throw std::runtime_error(std::string("Error, unknown token at ") + state);
				}
			}

			int lexpeek(parsectx & ctx, int eat = -1) {
				const char* preState = ctx.state;
				int rv = lexeat(ctx);
				if (rv == eat && rv != -1)
					return rv;
				ctx.state = preState;
				return rv;
			}


			std::shared_ptr<FTPL> parse_fun(parsectx & ctx, bool statement) {
				const char *srcstart = ctx.state;
				// TODO: add a skip-comment function
				while (*srcstart && isspace(*srcstart)) srcstart++;

				lexeat(ctx); // eat function token
				auto argscope = ctx.enter_scope(); // this enters the arg-scope.

				int argcount = 1;
				if (prep)
					argscope->names.emplace_back(rt->T_THIS, "this"); // dummy first arg to catch the "this"-arg

				std::string & tok = ctx.tok;
				auto id = std::make_pair(-1, std::string(""));
				// parse name part
				int name = lexpeek(ctx);
				if (name >= rt->USER) {
					lexeat(ctx);
					id = std::make_pair(name, ctx.tok);
				} else if (statement)
					throw std::runtime_error("expected function name but found " + tok);

				// next ensure we have parens for the argument list
				if (rt->LPAR != lexeat(ctx))
					throw std::runtime_error("expected ( after function name but found " + tok);
				// and parse the namelist until the end.
				if (rt->RPAR != lexpeek(ctx, rt->RPAR)) {
					while (true) {
						int arg = lexeat(ctx);
						if (arg<rt->USER)
							throw std::runtime_error("expected argument but found " + tok);
						if (prep)
							argscope->names.emplace_back(arg, tok);
						argcount++;
						int nxt = lexeat(ctx);
						if (rt->COMMA == nxt)
							continue;
						else if (rt->RPAR == nxt)
							break;
						else
							throw std::runtime_error("expected , or ) after argument but found " + tok);
					}
				}

				// store some previous state.
				GCTX *ogctx = ctx.gctx;
				// now setup things for our context.
				//printf(" ---- Begin fun....\n");
				GCTX gctx(rt);
				ctx.gctx = &gctx;
				argscope->maxstack = gctx.sp = 0;

				// check for a brace (the brace-pair will be parsed by the stmt parsing function)
				if (rt->LBRA != lexpeek(ctx))
					throw std::runtime_error("expected { after function arguments but found " + tok);
				// now parse the body
				stmt(ctx, -1);

				// always unconditinal return at end
				ctx.add(OPC::URETURN);

				// restore previous state
				ctx.gctx = ogctx;

				argscope->maxnames = argscope->names.size();

				ctx.leave_scope();
				//printf(" ---- End fun....\n");
				//auto outs = std::string(srcstart, ctx.state);
				//auto ptr = new FTPL(rt, id, argcount, argscope, gctx, std::move(outs));
				//return ctx.prep ? nullptr : ptr;
				return ctx.prep ? nullptr : std::make_shared<FTPL>(rt, id, argcount, argscope, gctx, std::string(srcstart, ctx.state));
			}
			std::unique_ptr<spres> expr(int & tt, parsectx & ctx, int prec) {
				std::string & tok = ctx.tok;
				double dnum;
				tt = lexeat(ctx, &dnum);
				if (tt == -1)
					return nullptr;
				auto stack = std::make_unique<spres>(ctx, 1);

				//OP cont = nullptr; OP val = nullptr; &cont, &val, 
				auto deref = [&]()->void {
					if (stack->sz != 1) {
						// TODO: do something here?!
						abort();
					}
				};
				if (tt >= rt->USER) {
					auto info = ctx.access(tt);
					// These instructions will be different between passes.. they cannot be differently sized!
					if (info.first == -1) {
						size_t idx = std::distance(ctx.gctx->globals.begin(), std::find(ctx.gctx->globals.begin(), ctx.gctx->globals.end(), tt));
						if (idx == ctx.gctx->globals.size())
							ctx.gctx->globals.push_back(tt);
						//printf("GlobLoad: %d %d\n",idx,tt);
						ctx.add(OPC::LOADGLOBAL, stack->off, idx);
					} else {
						ctx.add(OPC::LOADSSYM, stack->off, (info.second) | (info.first << 10));
					}
				} else if (rt->STR == tt) {
					ctx.add(OPC::LOADLIT, stack->off, gctx->addlit(tok)); // TODO: do value index..
																		//ctx.add(box<double>(dnum));
				} else if (tt == rt->DNUM) {
					ctx.add(OPC::LOADDOUBLE, stack->off, 0, false, &dnum);
				} else if (tt == rt->LPAR) {
					stack.reset();
					stack = expr(tt, ctx, 0);
					tt = lexeat(ctx);
					if (tt != rt->RPAR)
						throw std::runtime_error("Unexpected token in parenthised expression, wanted ) to end it but found " + ctx.tok);
				} else {
					throw std::runtime_error("unknown primtok " + ctx.tok);
				}

			operatorloop: while (-1 != (tt = lexpeek(ctx))) {
				// don't handle operators with too low precendce when doing a higher level.
				if (tt < prec)
					break;
				const std::tuple<int, int, OPC> binops[] = {
					{ rt->ADD,rt->MUL,OPC::ADD },
					{ rt->T_SUB,rt->MUL,OPC::SUB },
					{ rt->MUL,rt->T_MOD + 1,OPC::MUL },
					{ rt->T_DIV,rt->T_MOD + 1,OPC::DIV },
					{ rt->T_MOD,rt->T_MOD + 1,OPC::MOD },
					{ rt->EQ,rt->LT,OPC::EQ },
					{ rt->NEQ,rt->LT,OPC::NEQ },
					{ rt->LT,rt->GT + 1,OPC::LT },
					{ rt->LEQ,rt->GT + 1,OPC::LEQ },
					{ rt->GEQ,rt->GT + 1,OPC::GEQ },
					{ rt->GT,rt->GT + 1,OPC::GT }
				};
				for (int boidx = (sizeof(binops) / sizeof(binops[0])) - 1;boidx != -1;boidx--)
					if (std::get<0>(binops[boidx]) == tt) {
						lexeat(ctx); // eat the operator
						deref(); // deref any references since we only want a value
						std::unique_ptr<parsectx::spres> sub = expr(tt, ctx, std::get<1>(binops[boidx]));
						if (sub) {
							assert(stack->off + 1 == sub->off);
							ctx.add(std::get<2>(binops[boidx]), stack->off);
						}
						goto operatorloop;
					}
				// not a simple binary operator...
				if (rt->LPAR == tt) {
					// eat ( first
					lexeat(ctx);
					std::unique_ptr<parsectx::spres> thstmp;
					std::vector<std::unique_ptr<parsectx::spres>> targs;
					if (stack->sz == 1) {
						// no context given, only a fun-value
						// reserve one slot for the "this" value to be loaded as null
						thstmp = std::make_unique<parsectx::spres>(ctx, 1);
						ctx.add(OPC::LOADNULL, thstmp->off);
					} else {
						abort();
					}
					if (rt->RPAR == (tt = lexpeek(ctx, rt->RPAR))) {
					} else {
						while (true) {
							auto sub = expr(tt, ctx, 0);
							targs.push_back(std::move(sub));
							tt = lexeat(ctx);
							if (tt == -1)
								return nullptr;
							else if (rt->COMMA == tt)
								continue;
							else if (rt->RPAR == tt)
								break;
							else
								throw std::runtime_error("Unexpected token in funcall, wanted ) or , but found " + ctx.tok);
						}
					}
					// stack reserved for fun-call
					ctx.add(OPC::INVOKE, stack->off, targs.size() + 1);
				} else {
					break;
				}
			}
						  deref();
						  return std::move(stack);
			}
			int stmt(parsectx & ctx, int fnLvl) {
				std::string & tok = ctx.tok;

				int tt = lexpeek(ctx);
				if (tt == -1)
					return -1;
				if (rt->LBRA == tt) {
					lexeat(ctx);
					while (true) {
						if (rt->RBRA == lexpeek(ctx, rt->RBRA))
							break;
						if (-1 == stmt(ctx, fnLvl + 1))
							throw std::runtime_error("Premature EOF");
					}
					return 1;
				} else if (rt->T_FUNCTION == tt && fnLvl == 0) {
					auto fn = parse_fun(ctx, true);
					gctx->addfn(std::move(fn));
					//if (ctx.curscope) {
					//	abort();
					//} else if (!ctx.prep) {
					//	// we're at the toplevel eval, just box the fun-value directly!
					//	VAL v = box(std::move(funinst(nullptr, fn)));
					//	// we assign separately since we don't want bogus ptr slots in the map
					//	global[fn->id.first] = v;
					//}
					// ?? change to prefixed mkfun/setglob?
					// other options? add a fn-list to the parent-fun?
					// ^- very practical, doesn't disturb code-gen.
					return 1;
				} else if (rt->T_RETURN == tt) {
					lexeat(ctx);
					std::unique_ptr<parsectx::spres> rexp;
					if (rt->SEMICOLON == lexpeek(ctx, rt->SEMICOLON)) {
						// no return expr, just generate a null value.
						rexp = std::make_unique<parsectx::spres>(ctx, 1);
						ctx.add(OPC::LOADNULL);
					} else {
						rexp = expr(tt, ctx, 0);
						if (rt->SEMICOLON != lexeat(ctx))
							throw std::runtime_error("expected ; but found " + tok);
					}
					ctx.add(OPC::RETURN, rexp->off);
					return 1;
				} else if (rt->T_IF == tt) {
					// we need 2 gotos... IF fail target (else or end), tcode-end-target (post else or directly after?)
					lexeat(ctx);
					void* if_fail_key = (void*)ctx.state; // use left paren token ptr as if_fail_key
					if (rt->LPAR != lexpeek(ctx))
						throw std::runtime_error("expected ( after if but found " + tok);
					auto condslot = expr(tt, ctx, 0);
					// right parenthesis will have been parsed as part of the previous expression!
					ctx.add(OPC::FGOTO, condslot->off, ctx.cgmem[if_fail_key], true);
					// data in slot consumed!
					condslot.reset();
					stmt(ctx, fnLvl + 1);
					if (rt->T_ELSE == lexpeek(ctx)) {
						void* fcode_end_key = (void*)ctx.state; // use left paren token ptr as tcode_end_key
						lexeat(ctx);
						ctx.add(OPC::GOTO, 0, ctx.cgmem[fcode_end_key], true);
						ctx.cgmem[if_fail_key] = gctx->label();
						int rbr = stmt(ctx, fnLvl + 1);
						ctx.cgmem[fcode_end_key] = gctx->label();
					} else {
						ctx.cgmem[if_fail_key] = gctx->label();
					}
					return 1;
				} else if (rt->T_WHILE == tt) {
					void* while_stop_key = (void*)ctx.state; // use while token ptr as while_stop_key
					lexeat(ctx);
					if (rt->LPAR != lexpeek(ctx))
						throw std::runtime_error("expected ( after while but found " + tok);
					int retrypos = gctx->label();
					auto cond = expr(tt, ctx, 0);
					// right parenthesis will have been parsed as part of the expression!
					ctx.add(OPC::FGOTO, cond->off, ctx.cgmem[while_stop_key], true);
					cond.reset();
					// do body
					// TODO: support continue/break
					stmt(ctx, fnLvl + 1);
					// and goto start
					ctx.add(OPC::GOTO, 0, retrypos, true);
					// we should land here if the while cond fails
					ctx.cgmem[while_stop_key] = gctx->label();
					return 1;
				}
				auto op = expr(tt, ctx, 0);
				if (rt->SEMICOLON != lexeat(ctx))
					throw std::runtime_error("expected ; but found " + tok);
				return tt;
			}

			std::shared_ptr<FTPL> parse_top(const std::string & script, const std::string & src) {
				// eval runs the parser twice to build a scope chain during the first pass.
				//printf("---- Start of pass 1----------\n");
				state = script.c_str();
				topscope = std::make_shared<scopeinfo>();
				topscope->top = true;
				gctx->sp = 0;
				// run first pass of the parser.
				while (-1 != stmt(*this, 0)) {}

				//printf("---- End of pass 1----------\n");
				assert(gctx->sp == 0);
				assert(!curscope);
				// reset some things before the second iteration of the parser.
				prep = false;           // turn off the prep-flag, this will cause memory to be allocated generated.
				gctx->reset();
				state = script.c_str(); // reset the parser ptr
				tok.clear();            // clear the parsing token
				curscope = nullptr;     // and ensure that we are at the toplevel scope (should not be unless there is parser bugs)
				// also before the second pass we pre-dirty all scopes that has non-local variable accesses.
				for (auto& acc : accesses) {
					access(acc.second, acc.first);
				}
				// now run second pass o the parser
				while (-1 != stmt(*this, 0)) {}
				add(OPC::URETURN);
				return std::make_shared<FTPL>(
					rt,
					std::make_pair(-1, std::string()),
					0,
					topscope,
					*gctx,
					std::string(script));
			}


		};

		
		template<class RC, class EC>
		inline void do_cmp(INTERNALVALUE *slots, int off, const RC& rcmp, const EC& ecmp) {
			valuebase *lp = to_ptr(slots[off+0]);
			valuebase *rp = to_ptr(slots[off+1]);
			std::string * ls = lp ? lp->get<std::string>() : nullptr, *rs = rp ? rp->get<std::string>() : nullptr;
			bool * lb = lp ? lp->get<bool>() : nullptr, *rb = rp ? rp->get<bool>() : nullptr;
			if (ls && rs) {
				slots[off] = rcmp(*ls, *rs) ? V_TRUE : V_FALSE;
			} else if (lb && rb) {
				// do a truthiness comparasion
				slots[off] = ecmp(truthy(slots[off + 0]), truthy(slots[off + 1])) ? V_TRUE : V_FALSE;
			} else if (lp && rp) {
				// both objects!!
				slots[off] = ecmp(lp, rp) ? V_TRUE : V_FALSE;
//			} else { // only do this for both being doubles?
//				double lnum = unbox<double>(slots[off+0]), rnum = unbox<double>(slots[off+1]);
//				slots[off] =rcmp(lnum, rnum)?V_TRUE:V_FALSE;
			} else {
				throw std::exception("Tried comparing a number and a non-number value!");
			}
		}


	public:
		template<class T>
		void set(const std::string& id, T&& o) {
			int iid = atok(id.c_str());
			INTERNALVALUE v = box(std::move(o));
			global[iid] = v;
		}
		void eval(const std::string & script, const std::string & src = "<EVAL>") {
			parsectx<runtime,fungenctx, funtpl> ctx(this);
			fungenctx gctx(this);
			ctx.gctx = &gctx;
			std::shared_ptr<funtpl> efn = ctx.parse_top(script,src);
			// with the code generated, make a dummy fun and invoke it!
			{
				// TODO: efin needs to be able to mark the literals...
				funinst efin(nullptr, efn);
				efin.invoke(*this,0, nullptr);
			}
		}

		runtime(size_t heapsize=1000000) {
			size_t semisize = heapsize / 2;
			for (int i = 0;i < 2;i++) {
				heap[i].ptr = new char[semisize];
				heap[i].sz = semisize;
				hptr[i] = heap[i].ptr - i; // hptr's are a bit offset to simplify 
			}
			V_NULL = box((void*)nullptr);
			V_TRUE = box(true);
			V_FALSE = box(false);
			toksort();
			// end of built in tokens, register ID of first user-tok
			USER = nextTok;
		}
	};

	inline INTERNALVALUE runtime::funinst::invoke_impl(runtime *rt, funinst *fi, int argc, INTERNALVALUE *args) {
		auto ftpl=fi->qp;
		// create local scope if the functions scope is local.
		bool is_local = ftpl->info->local;
		INTERNALVALUE* localdata = is_local ? (INTERNALVALUE*)alloca(sizeof(INTERNALVALUE)*ftpl->info->maxnames) : nullptr;
		scope local(is_local?localdata:nullptr, ftpl->info.get());
		scope * cur = &local;
		if (!is_local) {
			abort(); // TODO, allocate a heap based scope instead! (hmm.... it could be feasible to do it just via a simple moveto since the logic is there...)
		}
		// TODO: varargs
		int copy_arg_count = std::min(argc, ftpl->argcount);
		for (int i = 0;i < copy_arg_count;i++) {
			cur->nslots[i] = args[i];
		}
		// install functions
		for (size_t i = 0;i < ftpl->functions.size();i++) {
			auto & fn = ftpl->functions[i];
			// TODO: scope
			INTERNALVALUE v = rt->box(funinst(nullptr, fn));
			if (ftpl->info->top) {
				// we assign separately since we don't want bogus ptr slots in the map
				rt->global[fn->id.first] = v;
			} else {
				abort();
			}
		}
		{
			auto ops = ftpl->code.data();
			int mstack = ftpl->info->maxstack;
			INTERNALVALUE* stack = (INTERNALVALUE*)alloca(sizeof(INTERNALVALUE)*mstack);
			frame cframe(rt, &cur,mstack, stack);
			INTERNALVALUE rv;
			rv.uval = 0;
			while (true) {
				int eop = *ops++;
				OPC op;
				int ARG0 = (eop >> 8) & 0xff;
				int ARG1 = (eop >> 16);
				switch (op = OPC(eop & 0xff)) {
				case OPC::URETURN: {
					rv = rt->V_NULL;
					goto eofun;
				}
				case OPC::RETURN: {
					rv = stack[ARG0];
					goto eofun;
				}
				case OPC::FGOTO: {
					if (!rt->truthy(stack[ARG0])) {
						ops += ARG1;
					}
					continue;
				}
				case OPC::GOTO: {
					ops += ARG1;
					continue;
				}
				case OPC::LOADNULL: {
					stack[ARG0] = rt->V_NULL;
					continue;
				}
				case OPC::INVOKE: {
					auto val = rt->invoke(stack[ARG0], ARG1, stack + ARG0 + 1);
					stack[ARG0] = val;
					continue;
				}
				case OPC::LOADDOUBLE: {
					stack[ARG0].uval = (((uint32_t)ops[0]) | (((int64_t)ops[1]) << 32));
					ops += 2;
					continue;
				}
				case OPC::LOADGLOBAL: {
					INTERNALVALUE* p = ftpl->globals[ARG1];
					stack[ARG0] = *p;
					continue;
				}
				case OPC::LOADSSYM: {
					scope* fs = cur;
					int up = ARG1>>10;
					while (up--) {
						fs = fs->parent;
					}
					stack[ARG0] = fs->nslots[ARG1&0x3ff];
					continue;
				}
				case OPC::LOADLIT:
					stack[ARG0] = ftpl->literals[ARG1].value;
					continue;
				case OPC::SUB:
					stack[ARG0] = rt->box(stack[ARG0].dval - stack[1 + ARG0].dval);
					continue;
				case OPC::MUL:
					stack[ARG0] = rt->box(stack[ARG0].dval * stack[1 + ARG0].dval);
					continue;
				case OPC::DIV:
					stack[ARG0] = rt->box(rt->unbox<double>(stack[ARG0]) / rt->unbox<double>(stack[1 + ARG0]));
					continue;
				case OPC::MOD:
					stack[ARG0] = rt->box(std::fmod(rt->unbox<double>(stack[ARG0]), rt->unbox<double>(stack[1 + ARG0])));
					continue;
#define NANOES__RUNTIME__RUNBLOCK__CMP(EOP,ROP) \
					if (stack[ARG0].uval<0xffffffff00000000LL&&stack[1 + ARG0].uval<0xffffffff00000000LL ) { \
						stack[ARG0]=( stack[ARG0].dval EOP stack[1 + ARG0].dval )?rt->V_TRUE:rt->V_FALSE; \
					} else { \
						rt->do_cmp(stack,ARG0,[](auto& a, auto& b) { return a EOP b; }, [](const auto& a, const auto& b)->bool { ROP }); \
					}
				case OPC::LT:
					NANOES__RUNTIME__RUNBLOCK__CMP(<, throw std::runtime_error("Cannot < compare a bool");)
					continue;
				case OPC::LEQ:
					NANOES__RUNTIME__RUNBLOCK__CMP(<= , throw std::runtime_error("Cannot <= compare a bool");)
					continue;
				case OPC::GEQ:
					NANOES__RUNTIME__RUNBLOCK__CMP(>= , throw std::runtime_error("Cannot >= compare a bool");)
					continue;
				case OPC::GT:
					NANOES__RUNTIME__RUNBLOCK__CMP(>, throw std::runtime_error("Cannot > compare a bool");)
					continue;
				case OPC::EQ:
					NANOES__RUNTIME__RUNBLOCK__CMP(== , return a == b;)
					continue;
				case OPC::NEQ:
					NANOES__RUNTIME__RUNBLOCK__CMP(!= , return a != b;)
					continue;
				case OPC::ADD: {
					if ( (stack[ARG0].uval| stack[ARG0+1].uval)<(uint64_t)0xffffffff00000000LL) {
						// fast-path if both numbers are far away from our NANTAG-pattern (should be majority of the time!)
						double dlv = stack[ARG0].dval, drv = stack[ARG0 + 1].dval;
						stack[ARG0] = rt->box(dlv + drv);
					} else {
						valuebase *lp = rt->to_ptr(stack[ARG0]);
						valuebase *rp = rt->to_ptr(stack[ARG0 + 1]);
						std::string * ls = lp ? lp->get<std::string>() : nullptr, *rs = rp ? rp->get<std::string>() : nullptr;
						if (ls || rs) {
							std::string sl = rt->unbox<std::string>(stack[ARG0]), sv = rt->unbox<std::string>(stack[ARG0 + 1]);
							std::string out = sl + sv;
							stack[ARG0] = rt->box<std::string>(std::move(out));
						} else {
							double dlv = stack[ARG0].dval, drv = stack[ARG0 + 1].dval;
							stack[ARG0] = rt->box(dlv + drv);
						}
					}
					continue;
				}
				case OPC::LOADSPROP :
					// obj from stack[ARG0] (need co-erce to obj before this call or number-proto-redir?)
					// if (!LOCALPTR)
					//   slowpath(&localptr);
					// else
					//   if (localptr[0]==hiddenclz) {
					//      stack[ARG0]=obj_as_cp+localptr[1];
					//   } else {
					//      for (i=2;;i+=2) {
					//         if (localptr[i]==hiddenclz) { swaplocalptr_i_and_0(); stack[ARG0]=obj_as_cp+localptr[i]; break;  }
					//         if (i<MAXLOCALPTR) continue;
					//         SLOWPATH(); break;
					//      }
					//      
					//   }
					//   INTRP( op+=LOCALPTRSIZE )
					throw std::runtime_error("Bad state " + std::to_string(int(op)));
				default:
					throw std::runtime_error("Bad state " + std::to_string(int(op)));
				}
			}
		eofun:
			return rv;
		}
	}

};


#endif // !INCLUDED_NANOES_HPP
