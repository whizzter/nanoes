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

	template<class T>
	intptr_t id_all() {
		static int id;
		return intptr_t(&id);
	}
	template<>
	intptr_t id_all<int>() {
		return 1;
	}
	template<>
	intptr_t id_all<double>() {
		return 2;
	}

	template<class T>
	intptr_t id() {
		return id_all<std::remove_const<std::remove_reference<T>::type>::type>();
	}

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
			int out = nextTok++;
			toks.emplace_back(out, v);
			return out;
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
				toksort();
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

		enum class STATUS {
			OK,
			RETURN,
			BREAK,
			CONTINUE
		};

		// Values are stored as NUN tagged 64 bit ints
		// If the high 32bits are clear then we have an object reference
		// where the lowest bit is a semi-space index and the rest is an offset into that space.
		// otherwise the 64bit number is an double that gets it's value by XOR:ing the specific NAN pattern.
		typedef int64_t VAL;
		static constexpr const int64_t NUNPAT = 0xfffff00f00000000LL;
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

		inline VAL invoke(runtime &rt,VAL fn, int arg, VAL* args) {
			if (valuebase *p = to_ptr(fn)) {
				return p->invoke(rt, arg, args);
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
			int max=0;
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
			int sp;

			struct scares {
				scope **mscope;
				int off;
				int sz;
				scares(scope**inscope, int x) : mscope(inscope),off((*inscope)->sp),sz(x) {
					(*mscope)->sp += sz;
				}
				~scares() {
					(*mscope)->sp -= sz;
				}
				VAL& operator[](int idx) {
#ifdef _DEBUG
					if (off + idx < 0 || off + idx >= (*mscope)->info->max)
						throw std::runtime_error("Error, out of bounds!");
#endif
					return (*mscope)->slots[off + idx];
				}
			};

			scope(VAL* islots, scopeinfo *in_info) : slots(islots), info(in_info),sp(in_info->names.size()) {
				if (islots)
					for (int i = 0;i < in_info->max;i++)
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

		struct opbase {
			virtual STATUS invoke(scope **,int roff) = 0;
			virtual int sym() { return -1; }
			virtual ~opbase() {}
			virtual int tsize() { return 0; }
		};
		typedef std::unique_ptr<opbase> OP;
		int tsize_block(std::vector<OP>& block) {
			return std::accumulate(block.begin(), block.end(), 0, [](auto prev, auto& or ) {
				return std::max(prev, or?or->tsize():0);
			});
		}
		int tsize_op(OP& op) {
			if (!op)
				return 0;
			return op->tsize();
		}

		VAL V_TRUE;  // TODO: set and collect
		VAL V_FALSE; // TODO: set and collect
		VAL V_NULL;  // TODO: set and collect

		template<class T>
		void touch(T* &v) {
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
			if (!heap[aHeap ^ 1].owns(old))
				throw std::exception("target not in old-gen");
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
		void touch(VAL& vr) {
			// only need to work on pointers, so change to a pointer and then rebase it..
			if (auto *p = to_ptr(vr)) {
				touch(p);
				vr=ptr_to_val(p);
			}
		}
		void collect(scope ** ss) {
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


		class nfun {
			friend runtime;
			runtime *rt;
			std::pair<int, std::string> id;
			std::shared_ptr<scopeinfo> info;
			int argcount;
			std::vector<OP> top;
		public:
			nfun(runtime *in_rt,std::pair<int, std::string> in_id, int in_argcount, std::shared_ptr<scopeinfo>& in_si, std::vector<OP> && in_top) : rt(in_rt),id(in_id), argcount(in_argcount), info(in_si), top(std::forward<std::vector<OP>>(in_top)) {}
			VAL invoke(int argc, VAL * args);
		};

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
			VAL invoke(runtime &rt, std::shared_ptr<nfun> &f, int argc, VAL *args) {
				return f->invoke(argc, args);
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
		inline typename std::enable_if< std::is_base_of<valuebase, T>::value, VAL>::type box(T && iv,scope** sa) {
			T *p = heap[aHeap].findfree<T>();
			if (!p) {
				collect(sa);
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
		inline typename std::enable_if< !std::is_base_of<valuebase, T>::value, VAL>::type box(T && iv, scope** sa) {
			return box(value<T>(std::move(iv)),sa);
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

		struct parsectx {
			const char *state;
			std::string tok;

			std::vector<std::pair<scopeinfo*, int>> accesses;
			std::vector<std::pair<const char*, std::shared_ptr<scopeinfo>>> scopes;
			std::shared_ptr<scopeinfo> curscope = nullptr;
			bool prep = true;

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

			template<class ROP>
			OP add_op(int tsize,ROP && op, int sym = -1) {
				if (prep) {
					return nullptr;
				} else {
					struct OPS : opbase {
						ROP data;
						int symbol;
						int mtsize;
						virtual int sym() { return symbol; }
						virtual int tsize() { return mtsize; }
						virtual STATUS invoke(scope ** vscope,int off) {
							return data(vscope,off);
						}
						OPS(int tsize, ROP&& ind, int insym) : mtsize(tsize), data(std::forward<ROP>(ind)), symbol(insym) {}
					};
					return std::make_unique<OPS>(tsize,std::forward<ROP>(op), sym);
				}
			}
		};


		static inline STATUS runblock(scope **scope,int roff, const std::vector<OP>& iops) {
			size_t sz = iops.size();
			auto ops = iops.data();
			for (size_t i = 0;i <sz;i++) {
				switch (ops[i]->invoke(scope,roff)) {
				case STATUS::OK:
					continue;
				case STATUS::RETURN:
					return STATUS::RETURN;
				case STATUS::BREAK:
					return STATUS::BREAK;
				case STATUS::CONTINUE:
					return STATUS::CONTINUE;
				default:
					abort();
				}
			}
			return STATUS::OK;
		}

		enum class OPC  {
			RETURN = 0,
			PUSHINT,   // one value (high 24 bits is the number)
			PUSHDOUBLE,  // 2 following values that are low/high bits
			PUSHLIT,     // high 24 bits is an literal stored inside the function!
			PUSHGLOB,    // ID as key? (hmm...)
			PUSHSCOPE,   // (12:12 bits) for up/offset
			INVOKE,      // high 24 bits is argcount
			ADD, SUB, MUL, DIV, MOD, // pops 2 and computes something.
			LT
		};
		VAL runblock(int32_t *ops) {
			scope **cur;
			while (true) {
				int eop = *ops++,op;
				switch ((OPC)(op = eop & 0xff)) {
				case OPC::PUSHDOUBLE :

				default:
					throw std::exception("Bad state!");
				}
			}
		}


		template<class MF>
		OP add_arith(parsectx&ctx, OP lop, OP rop, const MF & mf) {
			return ctx.add_op(2+std::max(tsize_op(lop),tsize_op(rop)),[this,lop = std::move(lop), rop = std::move(rop), mf](scope **scope,int roff) {
				scope::scares tmp(scope,2);
				lop->invoke(scope,tmp.off);
				rop->invoke(scope,tmp.off+1);
				double dlv = unbox<double>(tmp[0]), drv = unbox<double>(tmp[1]);
				(*scope)->slots[roff] = box(mf(dlv, drv));
				return STATUS::OK;
			});
		}
		template<class RC, class EC>
		OP add_cmp(parsectx&ctx, OP lop, OP rop, const RC & rcmp, const EC & ecmp) {
			return ctx.add_op(2 + std::max(tsize_op(lop), tsize_op(rop)),[this,lop = std::move(lop), rop = std::move(rop), rcmp, ecmp](scope ** scope,int roff) {
				scope::scares tmp(scope,2);
				lop->invoke(scope,tmp.off);
				rop->invoke(scope,tmp.off+1);

				valuebase *lp = to_ptr(tmp[0]);
				valuebase *rp = to_ptr(tmp[1]);
				std::string * ls = lp ? lp->get<std::string>() : nullptr, *rs = rp ? rp->get<std::string>() : nullptr;
				bool * lb = lp ? lp->get<bool>() : nullptr, *rb = rp ? rp->get<bool>() : nullptr;
				if (ls && rs) {
					(*scope)->slots[roff] = rcmp(*ls, *rs) ? V_TRUE : V_FALSE;
				} else if (lb && rb) {
					// do a truthiness comparasion
					(*scope)->slots[roff] =ecmp(truthy(tmp[0]),truthy(tmp[1]))?V_TRUE:V_FALSE;
				} else {
					double lnum = unbox<double>(tmp[0]), rnum = unbox<double>(tmp[1]);
					(*scope)->slots[roff] =rcmp(lnum, rnum)?V_TRUE:V_FALSE;
				}
				return STATUS::OK;
			});
		}



		std::shared_ptr<nfun> parse_fun(parsectx & ctx, bool statement) {

			auto argscope = ctx.enter_scope(); // this enters the arg-scope.

			int argcount = 1;
			argscope->names.emplace_back(T_THIS, "this"); // dummy first arg to catch the "this"-arg

			std::string & tok = ctx.tok;
			auto id = std::make_pair(-1, std::string(""));
			std::vector<OP> body;
			// parse name part
			int name = lexpeek(ctx);
			if (name >= USER) {
				lexeat(ctx);
				id = std::make_pair(name, ctx.tok);
			} else if (statement)
				throw std::runtime_error("expected function name but found " + tok);
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
			// check for a brace (the brace-pair will be parsed by the stmt parsing function)
			if (LBRA != lexpeek(ctx))
				throw std::runtime_error("expected { after function arguments but found " + tok);
			stmt(body, ctx, -1);
			ctx.leave_scope();
			argscope->max = 1+tsize_block(body) + argscope->names.size(); // reserve space for variables and tmps
			return std::make_shared<nfun>(this,id, argcount, argscope, std::move(body));
		}
		OP expr(int & tt, parsectx & ctx, int prec) {
			std::string & tok = ctx.tok;
			double dnum;
			tt = lexeat(ctx, &dnum);
			if (tt == -1)
				return nullptr;
			OP cont = nullptr;
			OP val = nullptr;
			auto deref = [this, &cont, &val, &ctx]()->void {
				if (cont) {
					val = ctx.add_op(tsize_op(cont)+tsize_op(val),[cont = std::move(cont), val = std::move(val)](scope **scope,int roff) {
						// TODO: deref from cont/val pairs
						abort();
						return STATUS::OK;
					});
					cont = nullptr;
				}
			};
			if (tt >= USER) {
				auto info = ctx.access(tt);
				if (info.first == -1) {
					if (global.find(tt) == global.end())
						global[tt] = NUNPAT;
					auto slotptr = &(global[tt]);
					val = ctx.add_op(0,[this, tt, slotptr](scope ** scope,int roff) {
						(*scope)->slots[roff] =(*slotptr);
						return STATUS::OK;
					}, tt);
				} else {
					val = ctx.add_op(0,[info](scope ** inscope,int roff) {
						int up = info.first;
						scope* scope = *inscope;
						while (up--) {
							scope = scope->parent;
						}
						(*inscope)->slots[roff] =scope->slots[info.second];
						return STATUS::OK;
					}, tt);
				}
			} else if (STR == tt) {
				val = ctx.add_op(0,[this,tok](scope ** scope,int roff) {
					(*scope)->slots[roff] = box(std::string(tok),scope);
					return STATUS::OK;
				});
			} else if (tt == DNUM) {
				val = ctx.add_op(0,[this,dnum](scope ** scope,int roff) {
					(*scope)->slots[roff] = box(dnum);
					return STATUS::OK;
				});
			} else if (tt == LPAR) {
				val = expr(tt, ctx, 0);
				tt = lexeat(ctx);
				if (tt != RPAR)
					throw std::runtime_error("Unexpected token in parenthised expression, wanted ) to end it but found " + ctx.tok);
			} else {
				throw std::runtime_error("unknown primtok " + ctx.tok);
			}

			while (-1 != (tt = lexpeek(ctx))) {
				// don't handle operators with too low precendce when doing a higher level.
				if (tt < prec)
					break;
				if (ADD == tt) {
					lexeat(ctx);
					deref(); // deref any references since we only want a value
					OP sub = expr(tt, ctx, MUL); // TODO: change precence req
					val = ctx.add_op(2+std::max(tsize_op(val),tsize_op(sub)),[this,val = std::move(val), sub = std::move(sub)](scope ** scope,int roff) {
						scope::scares tmp(scope, 2);
						val->invoke(scope,tmp.off);
						sub->invoke(scope,tmp.off+1);
						valuebase *lp = to_ptr(tmp[0]);
						valuebase *rp = to_ptr(tmp[1]);
						std::string * ls = lp?lp->get<std::string>():nullptr, *rs = rp?rp->get<std::string>():nullptr;
						if (ls || rs) {
							std::string sl = unbox<std::string>(tmp[0]), sv = unbox<std::string>(tmp[1]);
							std::string out = sl + sv;
							(*scope)->slots[roff]= box<std::string>(std::move(out),scope);
						} else {
							double dlv = unbox<double>(tmp[0]), drv = unbox<double>(tmp[1]);
							(*scope)->slots[roff] = box(dlv + drv);
						}
						return STATUS::OK;
					});
				} else if (T_DIV == tt || T_SUB == tt || T_MOD == tt || MUL == tt) {
					int op = tt;
					lexeat(ctx);
					deref();
					OP right = expr(tt, ctx, T_SUB==op?MUL:T_MOD+1);
					if (T_DIV == op) {
						val = add_arith(ctx, std::move(val), std::move(right), [](auto& a, auto& b) { return a / b; });
					} else if (T_SUB == op) {
						val = add_arith(ctx, std::move(val), std::move(right), [](auto& a, auto& b) { return a - b; });
					} else if (T_MOD == op) {
						val = add_arith(ctx, std::move(val), std::move(right), [](auto& a, auto& b) { return std::fmod(a, b); });
					} else if (MUL == op) {
						val = add_arith(ctx, std::move(val), std::move(right), [](auto& a, auto& b) { return a * b; });
					} else abort();
				} else if (LT == tt || LEQ == tt || GEQ == tt || GT == tt || EQ == tt || NEQ == tt) {
					int op = tt;
					lexeat(ctx);
					deref();
					OP right = expr(tt, ctx, (EQ==op || NEQ==op)?LT:GT+1);
					if (LT == op)
						val = add_cmp(ctx, std::move(val), std::move(right), [](auto& a, auto& b) { return a < b; }, [](const auto& a, const auto& b)->bool { throw std::runtime_error("Cannot < compare a bool"); });
					else if (LEQ == op)
						val = add_cmp(ctx, std::move(val), std::move(right), [](auto& a, auto& b) { return a <= b; }, [](const auto& a, const auto& b)->bool { throw std::runtime_error("Cannot <= compare a bool"); });
					else if (GEQ == op)
						val = add_cmp(ctx, std::move(val), std::move(right), [](auto& a, auto& b) { return a >= b; }, [](const auto& a, const auto& b)->bool { throw std::runtime_error("Cannot >= compare a bool"); });
					else if (GT == op)
						val = add_cmp(ctx, std::move(val), std::move(right), [](auto& a, auto& b) { return a > b; }, [](const auto& a, const auto& b)->bool { throw std::runtime_error("Cannot > compare a bool"); });
					else if (EQ == op)
						val = add_cmp(ctx, std::move(val), std::move(right), [](auto& a, auto& b) { return a == b; }, [](const auto& a, const auto& b) { return a == b; });
					else if (NEQ == op)
						val = add_cmp(ctx, std::move(val), std::move(right), [](auto& a, auto& b) { return a != b; }, [](const auto& a, const auto& b) { return a != b; });
					else abort();
				} else if (LPAR == tt) {
					// eat ( first
					lexeat(ctx);
					std::vector<OP> args;
					args.push_back(std::move(cont));
					args.push_back(std::move(val));
					cont = nullptr;
					if (RPAR == (tt = lexpeek(ctx, RPAR))) {
					} else {
						while (true) {
							OP sub = expr(tt, ctx, 0);
							args.push_back(std::move(sub));
							tt = lexeat(ctx);
							if (tt == -1)
								return nullptr;
							else if (COMMA == tt)
								continue;
							else if (RPAR == tt)
								break;
							else throw std::runtime_error("Unexpected token in funcall, wanted ) or , but found " + ctx.tok);
						}
					}
					int tsize = args.size() + tsize_block(args);
					//std::accumulate(args.begin(), args.end(), 0, [](auto& acc, auto&item) {return acc + tsize_op(item);});
					val = ctx.add_op(tsize,[this,args = std::move(args)](scope **scope,int roff) {
						scope::scares tmp(scope,args.size());

						if (args[0]) {
							abort(); // TODO: "this-invoke"
						} else {
							tmp[1] = V_NULL;
							args[1]->invoke(scope,tmp.off); // fetch fun.
						}
						for (size_t i = 2;i < args.size();i++) {
							args[i ]->invoke(scope,tmp.off+i);
						}

						(*scope)->slots[roff]=invoke(*this,tmp[0],args.size()-1, &tmp[1]);
						// scares destroys funvalue and args
						return STATUS::OK;
					});
				} else {
					break;
				}
			}
			deref();
			return std::move(val);
		}
		int stmt(std::vector<OP> & seq, parsectx & ctx, int fnLvl) {
			std::string & tok = ctx.tok;

			int tt = lexpeek(ctx);
			if (tt == -1)
				return -1;
			if (LBRA == tt) {
				lexeat(ctx);
				std::vector<OP> subblock;
				while (true) {
					if (RBRA == lexpeek(ctx, RBRA))
						break;
					if (-1 == stmt(subblock, ctx, fnLvl + 1))
						throw std::runtime_error("Premature EOF");
				}
				// TODO: check for scopes before just appending to the parent!
				for (size_t i = 0;i < subblock.size();i++)
					seq.push_back(std::move(subblock[i]));
				return 1;
			} else if (T_FUNCTION == tt && fnLvl == 0) {
				lexeat(ctx); // eat function token
				auto fn = parse_fun(ctx, true);
				if (ctx.curscope) {
					abort();
				} else {
					if (global.find(fn->id.first) == global.end())
						global[fn->id.first] = NUNPAT;
					auto slotptr = &(global[fn->id.first]);
					seq.insert(seq.begin(), ctx.add_op(0,[this, fn, slotptr](scope ** scope,int roff) {
						(*slotptr) = box(std::shared_ptr<nfun>(fn), scope);
						return STATUS::OK;
					}));
				}
				return 1;
			} else if (T_RETURN == tt) {
				lexeat(ctx);
				OP rexp;
				if (SEMICOLON == lexpeek(ctx, SEMICOLON)) {
				} else {
					rexp = expr(tt, ctx, 0);
					if (SEMICOLON != lexeat(ctx))
						throw std::runtime_error("expected ; but found " + tok);
				}
				seq.push_back(ctx.add_op(rexp?tsize_op(rexp):0,[this,rexp = std::move(rexp)](scope ** scope,int roff){
					if (rexp) {
						rexp->invoke(scope,roff);
					} else {
						(*scope)->slots[roff] = V_NULL;
					}
					return STATUS::RETURN;
				}));
				return 1;
			} else if (T_IF == tt) {
				lexeat(ctx);
				if (LPAR != lexpeek(ctx))
					throw std::runtime_error("expected ( after if but found " + tok);
				auto cond = expr(tt, ctx, 0);
				// right parenthesis will have been parsed as part of the expression!
				std::vector<OP> tcode;
				std::vector<OP> fcode;
				stmt(tcode, ctx, fnLvl + 1);
				if (T_ELSE == lexpeek(ctx, T_ELSE)) {
					int rbr = stmt(fcode, ctx, fnLvl + 1);
				}
				int tsize = 1 + std::max(tsize_op(cond), std::max(tsize_block(tcode), tsize_block(fcode)));
				seq.push_back(ctx.add_op(tsize, [this, cond = std::move(cond), tcode = std::move(tcode), fcode = std::move(fcode)](scope **scope, int roff) {
					scope::scares tmp(scope, 1);
					cond->invoke(scope, tmp.off);
					if (truthy(tmp[0])) {
						return runblock(scope, roff, tcode);
					} else {
						return runblock(scope, roff, fcode);
					}
				}));
				return 1;
			} else if (T_WHILE == tt) {
				lexeat(ctx);
				if (LPAR != lexpeek(ctx))
					throw std::runtime_error("expected ( after while but found " + tok);
				auto cond = expr(tt, ctx, 0);
				// right parenthesis will have been parsed as part of the expression!
				std::vector<OP> wcode;
				stmt(wcode, ctx, fnLvl + 1);
				int tsize = 1 + std::max(tsize_op(cond), tsize_block(wcode));
				seq.push_back(ctx.add_op(tsize, [this, cond = std::move(cond), wcode = std::move(wcode)](scope **scope, int roff) {
					scope::scares tmp(scope, 1);
					while (true) {
						cond->invoke(scope, tmp.off);
						if (!truthy(tmp[0]))
							break;
						STATUS code = runblock(scope, roff, wcode);
						// TODO: labels
						if (code == STATUS::OK || code == STATUS::CONTINUE)
							continue;
						else if (code == STATUS::BREAK)
							break;
						else if (code == STATUS::RETURN)
							return STATUS::RETURN;
						else abort();
					}
					return STATUS::OK;
				}));
				return 1;
			}
			auto op = expr(tt, ctx, 0);
			if (SEMICOLON != lexeat(ctx))
				throw std::runtime_error("expected ; but found " + tok);

			seq.push_back(std::move(op));
			return tt;
		}
	public:
		template<class T>
		void set(const std::string& id, T&& o) {
			while (true) {
				auto where = std::find_if(toks.begin(), toks.end(), [&](std::pair<int, std::string>& k) { return id == k.second; });
				if (where == toks.end()) {
					atok(id.c_str());
					toksort();
				} else {
					scope *s = nullptr;
					VAL v = box(std::move(o), &s);
					global[where->first] = v;
					return;
				}
			}
		}
		void eval(const std::string & script, const char *src = "<EVAL>") {
			// eval runs the parser twice to build a scope chain during the first pass.
			parsectx ctx;
			ctx.state = script.c_str();
			std::vector<OP> top;
			// run first pass of the parser.
			while (-1 != stmt(top, ctx, 0)) {}

			// reset some things before the second iteration of the parser.
			ctx.prep = false;           // turn off the prep-flag, this will cause memory to be allocated generated.
			top.clear();                // reset the top code vector.
			ctx.state = script.c_str(); // reset the parser ptr
			ctx.tok.clear();            // clear the parsing token
			ctx.curscope = nullptr;     // and ensure that we are at the toplevel scope (should not be unless there is parser bugs)
			// also before the second pass we pre-dirty all scopes that has non-local variable accesses.
			for (auto& acc : ctx.accesses) {
				ctx.access(acc.second, acc.first);
			}
			// now run second pass o the parser
			while (-1 != stmt(top, ctx, 0)) {}

			// with the code generated, make a dummy fun and invoke it!
			{
				auto si = std::make_shared<scopeinfo>();
				si->max = 1+tsize_block(top);
				nfun efn(this,std::make_pair(-1, std::string()), 0, si, std::move(top));
				efn.invoke(0, nullptr);
			}
		}

		runtime(size_t heapsize=1000000) {
			size_t semisize = heapsize / 2;
			for (int i = 0;i < 2;i++) {
				heap[i].ptr = new char[semisize];
				heap[i].sz = semisize;
				hptr[i] = heap[i].ptr - i; // hptr's are a bit offset to simplify 
			}
			scope *nscope = nullptr;
			V_NULL = box((void*)nullptr,&nscope);
			V_TRUE = box(true,&nscope);
			V_FALSE = box(false, &nscope);
			toksort();
			// end of built in tokens, register ID of first user-tok
			USER = nextTok;
		}
	};


	inline runtime::VAL runtime::nfun::invoke(int argc, VAL *args)
	{
		// create local scope if the functions scope is local.
		bool is_local = info->local;
		VAL* localdata = is_local ? (VAL*)alloca(sizeof(VAL)*info->max) : nullptr;
		scope local(is_local?localdata:nullptr, info.get());
		scope * fnscope = &local;
		if (!is_local) {
			abort(); // TODO, allocate a heap based scope instead! (hmm.... it could be feasible to do it just via a simple moveto since the logic is there...)
		}
		// TODO: varargs
		int copyargs = std::min(argc, argcount);
		for (int i = 0;i < copyargs;i++) {
			fnscope->slots[i] = args[i];
		}
		frame cur(rt, &fnscope);
		scope::scares tmp(&fnscope,1);
		runblock(&fnscope,tmp.off, top);
		return tmp[0];
	}

};


#endif // !INCLUDED_NANOES_HPP
