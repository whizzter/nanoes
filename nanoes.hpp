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

// direct-TDOP-compiler.
// problem 1: Symbols not visible? , enforce let or forbid var? problem when we need to late-discover!
//            trick... we could have a rewind position and just re-parse the entire code until the newly discovered snippet?!
//            ^- seems brilliant, could be a tad slow, could just make it 2-pass for function blocks (or per-script 2 pass with an initial pass of fun/binding recording)

namespace nanoes {
	class runtime;

	template<class T>
	static intptr_t id_all() {
		static int id;
		return intptr_t(&id);
	}
	template<>
	static intptr_t id_all<int>() {
		return 1;
	}
	template<>
	static intptr_t id_all<double>() {
		return 2;
	}

	template<class T>
	static intptr_t id() {
		return id_all<std::remove_const<std::remove_reference<T>::type>::type>();
	}

	// Values are stored as NUN tagged 64 bit ints
	// If the high 32bits are clear then we have an object reference
	// where the lowest bit is a semi-space index and the rest is an offset into that space.
	// otherwise the 64bit number is an double that gets it's value by XOR:ing the specific NAN pattern.
	typedef int64_t VAL;
	static constexpr const int64_t NUNPAT = 0xfffff00f00000000LL;

	class nesvalue {
		friend class runtime;
		runtime *rt;
		nesvalue *prev,*next;
		VAL value;
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
		nesvalue(runtime *in_rt) : prev(this), next(this), rt(in_rt), value(NUNPAT) {}
	public:
		nesvalue(const nesvalue& other) : rt(other.rt),value(other.value) {
			link(const_cast<nesvalue&>(other));
		}
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
		class gcbase;
		class valuebase;

		struct parsectx;
		struct opbase;

		static inline void* alignafter(size_t align, void *ip) {
			// find the next possible free position
			uintptr_t p = (uintptr_t)ip;
			// align our pointer properly.
			p = ((p + align - 1) / align)*align;
			return (void*)p;
		}

		// represent semi-spaces with the space structure.
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

		// this function indicates if the character is a valid symbol character
		int sym(int x, bool first = true) { return x == '_' || (first ? isalpha(x) : isalnum(x)); }

		// eat one token from the stream into the context
		int lexeat(parsectx & ctx, double *ddest = nullptr) {
			const char *& state = ctx.state;
			// skip over spaces (TODO: update line numbers?)
			while (std::isspace(*state)) state++;
			// search our sorted token+ID list by narrowing
			int rs = 0, re = toks.size() - 1;
			for (size_t sz = 0;rs <= re;) {
				// pre-token is smaller, go to next
				if (toks[rs].second.size() != sz && (toks[rs].second.size()<sz || toks[rs].second.c_str()[sz] < state[sz])) {
					rs++;
					continue;
				}
				// post-token is larger, go to next
				if (toks[re].second.size() != sz && (toks[re].second.size() < sz || toks[re].second.c_str()[sz] > state[sz])) {
					re--;
					continue;
				}
				// pre- and post- token are equal and we've found it. (also make sure we don't have an ID token and it continues!)
				if (rs == re && sz == toks[rs].second.size() && !(sym(state[0]) && sym(state[sz], false))) {
					ctx.tok.assign(state, sz);
					state += sz;
					return toks[rs].first;
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
				return STR;
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
				return DNUM;
			} else if (sym(*state)) {
				while (sym(*state, false))
					ctx.tok.push_back(*state++);
				int rv = atok(ctx.tok.c_str());
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

		// helper union for boxing/unboxing
		union VBOX {
			int64_t i64v;
			double dval;
		};

		class gcbase;
		//typedef std::function<void(VAL*, gcbase**)> toucher;
		typedef runtime toucher;

		std::map<int, VAL> global;

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
			// basic operations: prop-ref
			virtual VAL invoke(runtime &rt,int argc, VAL * args) = 0;
			virtual std::string to_string() = 0;
			virtual double to_num() = 0;
			virtual bool truthy() = 0;
		};

		inline valuebase* to_ptr(VAL v) {
			if (v & 0xffffffff00000000LL) {
				return 0;
			} else {
				int32_t i = (int32_t)v;
				return (valuebase*)(hptr[i & 1] + i);
			}
		}

		template<class T>
		inline T unbox(VAL v) {
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
		inline void* unbox<void*>(VAL v) {
			if (valuebase *p = to_ptr(v)) {
				return p;
			} else {
				return nullptr;
			}
		}
		template<>
		inline double unbox<double>(VAL v) {
			if (v & 0xffffffff00000000LL) {
				VBOX box;
				box.i64v = v^NUNPAT;
				return box.dval;
			} else {
				return to_ptr(v)->to_num();
			}
		}
		template<>
		inline std::string unbox<std::string>(VAL v) {
			if (valuebase *p = to_ptr(v)) {
				return p->to_string(); // TODO: do we want to move to GC strings?
			} else {
				return std::to_string(unbox<double>(v));
			}
		}
		template<>
		inline int unbox<int>(VAL v) {
			return (int)unbox<double>(v);
		}

		inline VAL invoke(VAL fn, int arg, VAL* args) {
			if (valuebase *p = to_ptr(fn)) {
				return p->invoke(*this, arg, args);
			} else {
				throw std::runtime_error(std::string("number ") + std::to_string(fn >> 1) + " is not a function");
			}
		}

		inline bool truthy(VAL v) {
			if (v == V_TRUE) return true;
			if (v == V_FALSE) return false;
			if (v == V_NULL) return false;
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
			std::shared_ptr<scopeinfo> parent = nullptr;
			std::vector<std::pair<int, std::string>> names;
			size_t max=0;
			bool local = true; // does variables in this scope escape? in that case we always need to allocate these scopes on the heap!
		};



		struct scope : public gcbase {
			virtual size_t size() {
				return (size_t)(((VAL*)alignafter(align(), (void*)sizeof(*this))) + info->max);
			}
			size_t runtime::scope::align()
			{
				return std::max(alignof(decltype(*this)), alignof(VAL));
			}
			virtual void moveto(void* dest) {
				scope * out = (scope*)dest;
				VAL* slots = (VAL*)alignafter(align(), out);
				memcpy(slots, this->slots, info->max);
				::new (out) scope(slots, info);
			}
			virtual void touch(toucher& to) {
				to.touch(parent);
				int max = info->max;
				for (int i = 0;i < max;i++)
					to.touch(slots[i]);
			}

			scope * parent = nullptr;

			scopeinfo *info;

			VAL* slots;

			scope(VAL* islots, scopeinfo *in_info) : slots(islots), info(in_info) {
				if (islots)
					for (size_t i = 0;i < in_info->max;i++)
						islots[i] = NUNPAT; // nunpat flips to zero
			}
		};
	
		struct frame;
		frame *activeframe= nullptr;
		struct frame {
			runtime *rt;
			frame *prev;
			scope** mscope;
			// TODO: move stack here?
			frame(runtime *in_rt, scope** in_scope) : mscope(in_scope), rt(in_rt), prev(rt->activeframe) {
				rt->activeframe = this;
			}
			~frame() {
				assert(rt->activeframe == this);
				rt->activeframe = prev;
			}
		};


		VAL V_TRUE;  // TODO: set and collect
		VAL V_FALSE; // TODO: set and collect
		VAL V_NULL;  // TODO: set and collect
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
		void touch(VAL& vr,bool allownew=false) {
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
			// mark the active stack frames.
			frame *t = activeframe;
			while (t) {
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
			VAL invoke_impl(runtime &rt, std::function<R(ARGS...)> & data, int argc, VAL * args, std::index_sequence<Is...>) {
				return rt.box(  data(   rt.unbox<ARGS>( args[Is] )...  ));
			}
			template<class R, class ... ARGS>
			VAL invoke(runtime &rt, std::function<R(ARGS...)> & data, int argc, VAL * args) {
				if (argc != sizeof...(ARGS)) {
					// TODO: location
					throw std::exception("Wrong script function arity");
				}
				return invoke_impl(rt,data, argc, args, std::make_index_sequence<sizeof...(ARGS)>{});
			}
			template<class T>
			VAL invoke(runtime &rt, T& ,int argc, VAL * args) {
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
			virtual VAL invoke(runtime &rt,int argc, VAL * args) {
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

		VAL ptr_to_val(void *p) {
			uintptr_t ip = (uintptr_t)p;
			auto dist = ip - ((uintptr_t)hptr[aHeap]);
			if (dist >= heap[aHeap].sz) {
				abort(); // TODO: fix so that the ptr_to_val macro can detect ptrs from the other heap (although that should never be done/used!)
			}
			return dist | aHeap;
		}


		template<class T>
		inline typename std::enable_if< std::is_base_of<valuebase, T>::value, VAL>::type box(T && iv) {
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
		inline typename std::enable_if< !std::is_base_of<valuebase, T>::value, VAL>::type box(T && iv) {
			return box(value<T>(std::move(iv)));
		}

		inline VAL box(double v) {
			VBOX vb;
			vb.dval = v;
			auto out = vb.i64v^NUNPAT;
			if (out& 0xffffffff00000000LL) {
				// if the pattern ok then return it
				return out;
			} else {
				// we had a NaN pattern that collided with our nunmask, return another nan
				return 0xfff8000000000000LL ^ NUNPAT;
			}
		}
		VAL box(int v) {
			return box((double)v);
		}

		enum class OPC : uint8_t {
			URETURN=0,
			RETURN,
			LOADINT,      // 2 value : next 12 is target-slot
			LOADDOUBLE,   // 3 values, 12 bits of target slot, 2 following values that are low/high bits
			LOADLIT,      // 1 value : 8bits op, 12bit target, ????
			LOADGLOB,     // 2 values :8bits op, 12bits target, 64bit ptr (only low 32bit loaded on 32bit arch)
			LOADNULL,     // 1 value  : 8bits op, 12bits target
			LOADTRUE,
			LOADFALSE,
			MOV,          // 2 values :8bits op, 12bits target, 12bits: next 12bits up (ignore masking initially!)
			INVOKE,       // 1 value : 8bit ops, 12bits base, 12 bits count
			ADD, SUB, MUL, DIV, MOD, // pops 2 and computes something. (8bit op, 12bit target, 12bit other)
			LT, LEQ, GEQ, GT, EQ, NEQ,
			FGOTO,        // 1 value : 8bit op, 12bit slot, 12bit jump offset
			GOTO          // 1 value : 8bit op, 12bit ????, 12bit jump offset
		};

		struct fungenctx {
			std::vector<int32_t> code;
			std::vector<nesvalue> literals;
			size_t sp = 0;
		};

		class funtpl {
			friend runtime;
			runtime *rt;
			std::pair<int, std::string> id;
			std::shared_ptr<scopeinfo> info;
			int argcount;
			std::vector<int32_t> code;
			std::vector<nesvalue> literals;
			std::string src;
		public:
			funtpl(runtime *in_rt, std::pair<int, std::string> in_id, int in_argcount, std::shared_ptr<scopeinfo>& in_si, fungenctx & in_gctx,std::string&& insrc)
				: rt(in_rt), id(in_id), argcount(in_argcount), info(in_si),code(std::move(in_gctx.code)),literals(in_gctx.literals),src(insrc)
			{}
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
			virtual bool truthy() {
				return true;
			}
			static inline VAL invoke_impl(runtime *rt, funinst *fi, int argc, VAL *args);
			virtual VAL invoke(runtime &rt, int argc, VAL *args) {
				return invoke_impl(&rt, this, argc, args);
			}

		};

		struct parsectx {
			runtime * rt;
			const char *state;
			std::string tok;

			parsectx(runtime *inrt):rt(inrt) {
//				assert(!rt->cureval);
//				rt->cureval = this;
			}
//			~parsectx() {
//				rt->cureval = nullptr;
//			}

			std::vector<std::pair<scopeinfo*, int>> accesses;
			std::vector<std::pair<const char*, std::shared_ptr<scopeinfo>>> scopes;
			std::shared_ptr<scopeinfo> topscope = nullptr;
			std::shared_ptr<scopeinfo> curscope = nullptr;
			bool prep = true;

			std::map<void*, int> cgmem;
			fungenctx * gctx;

			void add(int32_t v) {
				gctx->code.push_back(v);
			}
			void add(int64_t v) {
				add(int32_t(v));
				add(int32_t(v >> 32));
			}
			void addglob(const std::string& tok) {
				int id = rt->atok(tok.c_str());
				if (sizeof(int64_t) == sizeof(void*)) {
					int64_t addr = (int64_t)&(rt->global[id]);
					add(addr);
				} else if (sizeof(int32_t) == sizeof(void*)) {
					int32_t addr = (int32_t)&(rt->global[id]);
					add(addr);
				} else abort();
			}
			int addlit(const std::string& v) {
				int lid = gctx->literals.size();
				gctx->literals.emplace_back(rt->uroot);
				gctx->literals[lid].value = rt->box<std::string>(std::string(v));
				return lid;
			}
			void add(OPC a,int32_t b=0,int32_t c=0,bool reladdr=false) {
				if (reladdr) {
					c -= gctx->code.size() + 1;
				}
				add(((int)a) | (b << 8) | (c << 20));
			}
			int label() {
				return gctx->code.size();
			}

			struct spres {
				parsectx *ctx;
				size_t off;
				size_t sz;
				spres(parsectx &ictx, int x) : ctx(&ictx), off(ctx->gctx->sp), sz(x) {
					ctx->gctx->sp += sz;
					scopeinfo * ascope = ctx->curscope ? ctx->curscope.get() : ctx->topscope.get();
					ascope->max = std::max(ctx->gctx->sp, ascope->max);

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
						return std::make_pair(count, loc - from->names.begin());
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
		};

		
		template<class RC, class EC>
		inline void do_cmp(VAL *slots, int off, const RC& rcmp, const EC& ecmp) {
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

		std::shared_ptr<funtpl> parse_fun(parsectx & ctx, bool statement) {
			const char *srcstart = ctx.state;
			lexeat(ctx); // eat function token
			auto argscope = ctx.enter_scope(); // this enters the arg-scope.

			int argcount = 1;
			argscope->names.emplace_back(T_THIS, "this"); // dummy first arg to catch the "this"-arg

			std::string & tok = ctx.tok;
			auto id = std::make_pair(-1, std::string(""));
			// parse name part
			int name = lexpeek(ctx);
			if (name >= USER) {
				lexeat(ctx);
				id = std::make_pair(name, ctx.tok);
			} else if (statement)
				throw std::runtime_error("expected function name but found " + tok);

			// store a key for our namecount 
			void * namecountkey = (void*)ctx.state;
			// next ensure we have parens for the argument list
			if (LPAR != lexeat(ctx))
				throw std::runtime_error("expected ( after function name but found " + tok);
			// and parse the namelist until the end.
			if (RPAR != lexpeek(ctx, RPAR)) {
				while (true) {
					int arg = lexeat(ctx);
					if (arg<USER)
						throw std::runtime_error("expected argument but found " + tok);
					argscope->names.emplace_back(arg, tok);
					argcount++;
					int nxt = lexeat(ctx);
					if (COMMA == nxt)
						continue;
					else if (RPAR == nxt)
						break;
					else
						throw std::runtime_error("expected , or ) after argument but found " + tok);
				}
			}
			
			// store some previous state.
			fungenctx *ogctx= ctx.gctx;
			// now setup things for our context.
			fungenctx gctx;
			ctx.gctx = &gctx;
			gctx.sp = ctx.prep ? 0 : ctx.cgmem[namecountkey];
			argscope->max = gctx.sp;

			// check for a brace (the brace-pair will be parsed by the stmt parsing function)
			if (LBRA != lexpeek(ctx))
				throw std::runtime_error("expected { after function arguments but found " + tok);
			// now parse the body
			stmt(ctx, -1);

			// always unconditinal return at end
			ctx.add(OPC::URETURN);

			// restore previous state
			ctx.gctx = ogctx;

			// memorize the number of names so we can calculate accurate stack positions next time.
			if (ctx.prep)
				ctx.cgmem[namecountkey] = argscope->names.size();

			ctx.leave_scope();
			return std::make_shared<funtpl>(this,id, argcount, argscope, gctx,std::string(srcstart,ctx.state));
		}
		std::unique_ptr<parsectx::spres> expr(int & tt, parsectx & ctx, int prec) {
			std::string & tok = ctx.tok;
			double dnum;
			tt = lexeat(ctx, &dnum);
			if (tt == -1)
				return nullptr;
			auto stack = std::make_unique<parsectx::spres>(ctx, 1);
			
			//OP cont = nullptr; OP val = nullptr; &cont, &val, 
			auto deref = [&]()->void {
				if (stack->sz != 1) {
					// TODO: do something here?!
					abort();
				}
			};
			if (tt >= USER) {
				auto info = ctx.access(tt);
				if (info.first == -1) {
					ctx.add(OPC::LOADGLOB, stack->off);
					ctx.addglob(ctx.tok);
				} else {
					ctx.add(OPC::MOV,stack->off,info.second);
					ctx.add(info.first); // up stored
				}
			} else if (STR == tt) {
				ctx.add(OPC::LOADLIT, stack->off,ctx.addlit(tok)); // TODO: do value index..
				//ctx.add(box<double>(dnum));
			} else if (tt == DNUM) {
				ctx.add(OPC::LOADDOUBLE, stack->off);
				ctx.add(box(dnum));
			} else if (tt == LPAR) {
				stack.reset();
				stack=expr(tt, ctx, 0);
				tt = lexeat(ctx);
				if (tt != RPAR)
					throw std::runtime_error("Unexpected token in parenthised expression, wanted ) to end it but found " + ctx.tok);
			} else {
				throw std::runtime_error("unknown primtok " + ctx.tok);
			}

			operatorloop: while (-1 != (tt = lexpeek(ctx))) {
				// don't handle operators with too low precendce when doing a higher level.
				if (tt < prec)
					break;
				const std::tuple<int, int, OPC> binops[] = {
					{ ADD,MUL,OPC::ADD },
					{ T_SUB,MUL,OPC::SUB },
					{ MUL,T_MOD + 1,OPC::MUL },
					{ T_DIV,T_MOD + 1,OPC::DIV },
					{ T_MOD,T_MOD + 1,OPC::MOD },
					{ EQ,LT,OPC::EQ },
					{ NEQ,LT,OPC::NEQ },
					{ LT,GT + 1,OPC::LT },
					{ LEQ,GT + 1,OPC::LEQ },
					{ GEQ,GT + 1,OPC::GEQ },
					{ GT,GT + 1,OPC::GT }
				};
				for (int boidx = (sizeof(binops) / sizeof(binops[0])) - 1;boidx != -1;boidx--)
					if (std::get<0>(binops[boidx])==tt) {
						lexeat(ctx); // eat the operator
						deref(); // deref any references since we only want a value
						std::unique_ptr<parsectx::spres> sub = expr(tt, ctx, std::get<1>(binops[boidx]) );
						if (sub) {
							assert(stack->off+1 == sub->off);
							ctx.add(std::get<2>(binops[boidx]), stack->off);
						}
						goto operatorloop;
					}
				// not a simple binary operator...
				if (LPAR == tt) {
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
					if (RPAR == (tt = lexpeek(ctx, RPAR))) {
					} else {
						while (true) {
							auto sub = expr(tt, ctx, 0);
							targs.push_back(std::move(sub));
							tt = lexeat(ctx);
							if (tt == -1)
								return nullptr;
							else if (COMMA == tt)
								continue;
							else if (RPAR == tt)
								break;
							else
								throw std::runtime_error("Unexpected token in funcall, wanted ) or , but found " + ctx.tok);
						}
					}
					// stack reserved for fun-call
					ctx.add(OPC::INVOKE, stack->off, targs.size()+1);
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
			if (LBRA == tt) {
				lexeat(ctx);
				while (true) {
					if (RBRA == lexpeek(ctx, RBRA))
						break;
					if (-1 == stmt(ctx, fnLvl + 1))
						throw std::runtime_error("Premature EOF");
				}
				return 1;
			} else if (T_FUNCTION == tt && fnLvl == 0) {
				auto fn = parse_fun(ctx, true);
				if (ctx.curscope) {
					abort();
				} else {
					// we're at the toplevel eval, just box the fun-value directly!
					VAL v = box(std::move(funinst(nullptr, fn)));
					// we assign separately since we don't want bogus ptr slots in the map
					global[fn->id.first] = v;
				}
				return 1;
			} else if (T_RETURN == tt) {
				lexeat(ctx);
				std::unique_ptr<parsectx::spres> rexp;
				if (SEMICOLON == lexpeek(ctx, SEMICOLON)) {
					// no return expr, just generate a null value.
					rexp = std::make_unique<parsectx::spres>(ctx, 1);
					ctx.add(OPC::LOADNULL);
				} else {
					rexp = expr(tt, ctx, 0);
					if (SEMICOLON != lexeat(ctx))
						throw std::runtime_error("expected ; but found " + tok);
				}
				ctx.add(OPC::RETURN,rexp->off);
				return 1;
			} else if (T_IF == tt) {
				// we need 2 gotos... IF fail target (else or end), tcode-end-target (post else or directly after?)
				void* if_fail_key = (void*)ctx.state; // use if token ptr as if_fail_key
				lexeat(ctx);
				if (LPAR != lexpeek(ctx))
					throw std::runtime_error("expected ( after if but found " + tok);
				auto condslot = expr(tt, ctx, 0);
				// right parenthesis will have been parsed as part of the previous expression!
				ctx.add(OPC::FGOTO, condslot->off, ctx.cgmem[if_fail_key], true);
				// data in slot consumed!
				condslot.reset();
				stmt(ctx, fnLvl + 1);
				if (T_ELSE == lexpeek(ctx, T_ELSE)) {
					void* fcode_end_key = (void*)ctx.state; // use left paren token ptr as tcode_end_key
					ctx.add(OPC::GOTO, 0, ctx.cgmem[fcode_end_key], true);
					ctx.cgmem[if_fail_key] = ctx.label();
					int rbr = stmt(ctx, fnLvl + 1);
					ctx.cgmem[fcode_end_key] = ctx.label();
				} else {
					ctx.cgmem[if_fail_key] = ctx.label();
				}
				return 1;
			} else if (T_WHILE == tt) {
				void* while_stop_key = (void*)ctx.state; // use while token ptr as while_stop_key
				lexeat(ctx);
				if (LPAR != lexpeek(ctx))
					throw std::runtime_error("expected ( after while but found " + tok);
				int retrypos = ctx.label();
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
				ctx.cgmem[while_stop_key] = ctx.label();
				return 1;
			}
			auto op = expr(tt, ctx, 0);
			if (SEMICOLON != lexeat(ctx))
				throw std::runtime_error("expected ; but found " + tok);
			return tt;
		}
	public:
		template<class T>
		void set(const std::string& id, T&& o) {
			int iid = atok(id.c_str());
			VAL v = box(std::move(o));
			global[iid] = v;
		}
		void eval(const std::string & script, const char *src = "<EVAL>") {
			std::shared_ptr<funtpl> efn;
			{
				// eval runs the parser twice to build a scope chain during the first pass.
				parsectx ctx(this);
				ctx.rt = this;
				ctx.state = script.c_str();
				ctx.topscope = std::make_shared<scopeinfo>();
				fungenctx gctx;
				ctx.gctx = &gctx;
				gctx.sp = 0;
				// run first pass of the parser.
				while (-1 != stmt(ctx, 0)) {}

				assert(gctx.sp == 0);
				assert(!ctx.curscope);
				// reset some things before the second iteration of the parser.
				ctx.prep = false;           // turn off the prep-flag, this will cause memory to be allocated generated.
				gctx.code.clear();                // reset the top code vector.
				gctx.literals.clear();       // remove the literals to re-gen them
				ctx.state = script.c_str(); // reset the parser ptr
				ctx.tok.clear();            // clear the parsing token
				ctx.curscope = nullptr;     // and ensure that we are at the toplevel scope (should not be unless there is parser bugs)
				// also before the second pass we pre-dirty all scopes that has non-local variable accesses.
				for (auto& acc : ctx.accesses) {
					ctx.access(acc.second, acc.first);
				}
				// now run second pass o the parser
				while (-1 != stmt(ctx, 0)) {}
				ctx.add(OPC::URETURN);
				efn = std::make_shared<funtpl>(this, std::make_pair(-1, std::string()), 0, ctx.topscope, gctx,std::string(script));
			}
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

	inline VAL runtime::funinst::invoke_impl(runtime *rt, funinst *fi, int argc, VAL *args) {
		auto ftpl=fi->qp;
		// create local scope if the functions scope is local.
		bool is_local = ftpl->info->local;
		VAL* localdata = is_local ? (VAL*)alloca(sizeof(VAL)*ftpl->info->max) : nullptr;
		scope local(is_local?localdata:nullptr, ftpl->info.get());
		scope * cur = &local;
		if (!is_local) {
			abort(); // TODO, allocate a heap based scope instead! (hmm.... it could be feasible to do it just via a simple moveto since the logic is there...)
		}
		// TODO: varargs
		int copyargs = std::min(argc, ftpl->argcount);
		for (int i = 0;i < copyargs;i++) {
			cur->slots[i] = args[i];
		}
		//	return rt.runblock(qp,fnscope, this->code->code.data());
		{
			auto ops = ftpl->code.data();
			frame cframe(rt, &cur);
			VAL rv = NUNPAT;
			while (true) {
				int eop = *ops++;
				OPC op;
				int ARG0 = (eop >> 8) & 0xfff;
				int ARG1 = (eop >> 20);
				switch (op = OPC(eop & 0xff)) {
				case OPC::URETURN: {
					rv = rt->V_NULL;
					goto eofun;
				}
				case OPC::RETURN: {
					rv = cur->slots[ARG0];
					goto eofun;
				}
				case OPC::FGOTO: {
					if (!rt->truthy(cur->slots[ARG0])) {
						ops += ARG1;
					}
					continue;
				}
				case OPC::GOTO: {
					ops += ARG1;
					continue;
				}
				case OPC::LOADNULL: {
					cur->slots[ARG0] = rt->V_NULL;
					continue;
				}
				case OPC::INVOKE: {
					auto val = rt->invoke(cur->slots[ARG0], ARG1, cur->slots + ARG0 + 1);
					cur->slots[ARG0] = val;
					continue;
				}
				case OPC::LOADDOUBLE: {
					cur->slots[ARG0] = (((uint32_t)ops[0]) | (((int64_t)ops[1]) << 32));
					ops += 2;
					continue;
				}
				case OPC::LOADGLOB: {
					if (sizeof(int64_t) == sizeof(void*)) {
						VAL* p = (VAL*)(((uint32_t)ops[0]) | (((int64_t)ops[1]) << 32));
						ops += 2;
						cur->slots[ARG0] = *p;
					} else if (sizeof(int32_t) == sizeof(void*)) {
						VAL* p = (VAL*)*ops++;
						cur->slots[ARG0] = *p;
					} else abort();
					continue;
				}
				case OPC::MOV: {
					scope* fs = cur;
					int up = *ops++;
					while (up--) {
						fs = fs->parent;
					}
					cur->slots[ARG0] = fs->slots[ARG1];
					continue;
				}
				case OPC::LOADLIT:
					cur->slots[ARG0] = ftpl->literals[ARG1].value;
					continue;
				case OPC::SUB:
					cur->slots[ARG0] = rt->box(rt->unbox<double>(cur->slots[ARG0]) - rt->unbox<double>(cur->slots[1 + ARG0]));
					continue;
				case OPC::MUL:
					cur->slots[ARG0] = rt->box(rt->unbox<double>(cur->slots[ARG0]) * rt->unbox<double>(cur->slots[1 + ARG0]));
					continue;
				case OPC::DIV:
					cur->slots[ARG0] = rt->box(rt->unbox<double>(cur->slots[ARG0]) / rt->unbox<double>(cur->slots[1 + ARG0]));
					continue;
				case OPC::MOD:
					cur->slots[ARG0] = rt->box(std::fmod(rt->unbox<double>(cur->slots[ARG0]), rt->unbox<double>(cur->slots[1 + ARG0])));
					continue;
#define NANOES__RUNTIME__RUNBLOCK__CMP(EOP,ROP) \
					if ((0xffffffff00000000LL&cur->slots[ARG0])&&(0xffffffff00000000LL&cur->slots[1 + ARG0])) { \
						cur->slots[ARG0]=( rt->unbox<double>(cur->slots[ARG0]) EOP rt->unbox<double>(cur->slots[1 + ARG0]) )?rt->V_TRUE:rt->V_FALSE; \
					} else { \
						rt->do_cmp(cur->slots,ARG0,[](auto& a, auto& b) { return a EOP b; }, [](const auto& a, const auto& b)->bool { ROP }); \
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
					if ((0xffffffff00000000LL & cur->slots[ARG0]) && (0xffffffff00000000LL & cur->slots[1 + ARG0])) {
						goto op_add_double;
					}
					valuebase *lp = rt->to_ptr(cur->slots[ARG0]);
					valuebase *rp = rt->to_ptr(cur->slots[ARG0 + 1]);
					std::string * ls = lp ? lp->get<std::string>() : nullptr, *rs = rp ? rp->get<std::string>() : nullptr;
					if (ls || rs) {
						std::string sl = rt->unbox<std::string>(cur->slots[ARG0]), sv = rt->unbox<std::string>(cur->slots[ARG0 + 1]);
						std::string out = sl + sv;
						cur->slots[ARG0] = rt->box<std::string>(std::move(out));
					} else {
					op_add_double:
						double dlv = rt->unbox<double>(cur->slots[ARG0]), drv = rt->unbox<double>(cur->slots[ARG0 + 1]);
						cur->slots[ARG0] = rt->box(dlv + drv);
					}
					continue;
				}
				default:
					throw std::runtime_error("Bad state " + std::to_string(int(op)));
				}
			}
		eofun:
			// after returning we clear out the stack-values so we don't hold extra references. (make it exception safe?)
			for (size_t i = cur->info->names.size();i < cur->info->max;i++) {
				cur->slots[i] = NUNPAT;
			}
			return rv;
		}
	}

};


#endif // !INCLUDED_NANOES_HPP
