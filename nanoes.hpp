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

#ifdef _MSC_VER
#include <malloc.h>
#else
#include <alloca.h>
#endif

// direct-TDOP-compiler.
// problem 1: Symbols not visible? , enforce let or forbid var? problem when we need to late-discover!
//            trick... we could have a rewind position and just re-parse the entire code until the newly discovered snippet?!
//            ^- seems brilliant, could be a tad slow, could just make it 2-pass for function blocks (or per-script 2 pass with an initial pass of fun/binding recording)

namespace nanoes {
	// opaque typelist
	template<class ... T>
	struct typeset;

	//using deftypes = typeset<int, double, std::string, std::vector<intptr_t>, std::map<intptr_t, intptr_t>>;
	template<class T> struct tslen { static const int sz = 0; };
	template<class H, class ...R> struct tslen<typeset<H, R...>> { static const int sz = tslen<typeset<R...>>::sz + 1; };

	template<int I, class TS> struct tsinfo {
		void idx() {}
	};
	template<int I, class H, class ...R> struct tsinfo<I, typeset<H, R...>> : tsinfo<I + 1, typeset<R...>> {
		using parent = tsinfo<I + 1, typeset<R...>>;
		using parent::idx;
		constexpr int idx(H*) const { return I; };
	};

	template<class T>
	intptr_t id_all() {
		static int id;
		return intptr_t(&id);
	}

	template<class T>
	intptr_t id() {
		return id_all< std::remove_const<std::remove_reference<T>::type>::type >();
	}

	//template<class TS=typename deftypes>
	class runtime {
		class valuebase;
		struct slot;
		struct parsectx;

		// lexer/parser info and a complete symbol table.
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
		const int LPAR = atok("("), RPAR = atok(")"), DOT = atok(".");
		const int T_IF = atok("if"), T_ELSE = atok("else"), T_WHILE = atok("while"), T_FUNCTION = atok("function"), T_RETURN = atok("return"), T_THIS = atok("this");
		const int ID = atok("   _ID_"), DNUM = atok("   _dnum_"), STR = atok("   __str__");
		int USER;

		// try to derive a SLOTSIZE that can be useful for everything needed in the runtime
		static constexpr int SLOTSIZE = std::max(sizeof(std::string) * 2, std::max(sizeof(std::vector<double>) * 2, std::max(sizeof(double) * 5, sizeof(void*) * 5)));


		std::map<int, std::unique_ptr<valuebase>> global;

		void toksort() {
			std::sort(toks.begin(), toks.end(), [](auto a, auto b) {  return a.second < b.second; });
		}

		int sym(int x, bool first = true) { return x == '_' || (first ? isalpha(x) : isalnum(x)); }

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
			ctx.tok.clear();
			// ok not an known operator or identifier, either a literal
			if (*state == '\"') {
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
						// TODO: charcter escape codes.
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
				// TODO: fractions and exponent
				return DNUM; //intish ? INUM : DNUM;
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

		//using STATUS = void;

		class valuebase {
		protected:
			virtual intptr_t ty() = 0;// { abort(); }
		public:
			virtual ~valuebase() {}
			virtual size_t size() = 0;
			virtual void moveto(slot& data) = 0;// { abort(); }
			virtual void copyto(slot& data) = 0;// { abort(); }

			virtual void * get(void**p = nullptr) { abort(); }

			template<class T> inline T* get(T* x = nullptr) {
				if (ty() == id<T>())
					return (T*)get();
				else
					return nullptr;
			}
			inline void* ref(void**p = nullptr) {
				return get();
			}
			template<class T> inline T& ref(T* p = nullptr) {
				if (ty() == id<T>())
					return *(T*)get();
				else
					throw std::exception("nonmatching script type value");
			}
			// basic operations: prop-ref
			virtual void invoke(slot&, int argc, slot * args) = 0; // { abort(); }
			virtual std::string to_string() = 0;
			virtual double to_num() = 0;
			virtual bool truthy() = 0;
		};

		struct slot {
			alignas(std::max(alignof(double), alignof(void*))) char data[SLOTSIZE];

			template<class T, class ... ARG>
			inline void set(ARG&&... args);
			inline valuebase* val() { return (valuebase*)data; }
			inline slot();
			inline ~slot();
		};

		// 1st pass, create scopes and note scopes. Also register usages.
		// mid-pass, find variable targets and mark scopes non-local.
		// 2nd pass, just use scope-info to 
		struct scopeinfo {
			std::shared_ptr<scopeinfo> parent = nullptr;
			std::vector<std::pair<int, std::string>> names;
			bool local = true; // does variables in this scope escape? in that case we always need to allocate these scopes on the heap!
		};


		struct scope {
			scope * parent = nullptr;

			scopeinfo *info;
			int count; // some duplication here

			slot* slots;

			scope(slot* islots, int icount, scopeinfo *in_info) : slots(islots), count(icount), info(in_info) {}
		};

		struct opbase {
			virtual STATUS invoke(slot& res, scope *) = 0;
			virtual int sym() { return -1; }
			virtual ~opbase() {}
		};
		typedef std::unique_ptr<opbase> OP;


		class nfun {
			friend runtime;
			std::pair<int, std::string> id;
			std::shared_ptr<scopeinfo> info;
			int argcount;
			std::vector<OP> top;
		public:
			nfun(std::pair<int, std::string> in_id, int in_argcount, std::shared_ptr<scopeinfo>& in_si, std::vector<OP> && in_top) : id(in_id), argcount(in_argcount), info(in_si), top(std::forward<std::vector<OP>>(in_top)) {}
			void invoke(slot& out, int argc, slot * args);
		};

		template<class T>
		class value : public valuebase {
			T data;

			template<class R, class ... ARGS, size_t... Is>
			void invoke_impl(std::function<R(ARGS...)> & data, slot& dest, int argc, slot * args, std::index_sequence<Is...>) {
				auto rv = data(args[Is].val()->ref((ARGS*)nullptr)...);
				dest.set<R>(std::move(rv));
			}
			template<class R, class ... ARGS>
			void invoke(std::shared_ptr<std::function<R(ARGS...)>> & data, slot& dest, int argc, slot * args) {
				if (argc != sizeof...(ARGS)) {
					// TODO: location
					throw std::exception("Wrong script function arity");
				}
				invoke_impl(*data, dest, argc, args, std::make_index_sequence<sizeof...(ARGS)>{});
			}
			void invoke(std::shared_ptr<nfun> &f, slot &dest, int argc, slot *args) {
				f->invoke(dest, argc, args);
			}
			template<class T>
			void invoke(T & data, slot& dest, int argc, slot * args) {
				std::cout << "Invoking for ANY\n";
				dest.set<bool>(false);
			}

			double to_num(std::string & s) {
				throw std::runtime_error("Not a number : [" + s + "]");
			}
			double to_num(double d) {
				return d;
			}
			double to_num(int v) {
				return v;
			}
			template<class T>
			double to_num(T & v) {
				throw std::runtime_error(typeid(v).name() + std::string(" is not a number"));
			}

			bool truthy(int v) {
				return v != 0;
			}
			bool truthy(double v) {
				return v != 0;
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
			std::string to_string(int v) {
				return std::to_string(v);
			}
			std::string to_string(double v) {
				return std::to_string(v);
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
			virtual void *get(void **vp = nullptr) { return &data; }
			virtual void moveto(slot& dest) { dest.set<T>(std::move(*this)); } //new(dest) value(std::move(*this)); }
			virtual void copyto(slot& dest) { dest.set<T>(*this); } //new(dest) value(*this); }
			virtual void invoke(slot& dest, int argc, slot * args) {
				invoke(data, dest, argc, args);
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
			OP add_op(ROP && op, int sym = -1) {
				if (prep) {
					return nullptr;
				} else {
					struct OPS : opbase {
						ROP data;
						int symbol;
						virtual int sym() { return symbol; }
						OPS(ROP&& ind, int insym) : data(std::forward<ROP>(ind)), symbol(insym) {}
						virtual STATUS invoke(slot& res, scope * vscope) {
							return data(res, vscope);
						}
					};
					return std::make_unique<OPS>(std::forward<ROP>(op), sym);
				}
			}
		};


		static inline STATUS runblock(slot &out, scope *scope, const std::vector<OP>& iops) {
			size_t sz = iops.size();
			auto ops = iops.data();
			for (size_t i = 0;i <sz;i++) {
				switch (ops[i]->invoke(out, scope)) {
				case STATUS::OK:
					continue;
				case STATUS::RETURN:
					return STATUS::RETURN;
				case STATUS::BREAK: // todo
				case STATUS::CONTINUE: // todo
				default:
					abort();
				}
			}
			return STATUS::OK;
		}

		struct sscope : scope {
			sscope(slot* islots, int icount, scopeinfo *in_info = nullptr) : scope(islots, icount, in_info) {
				for (int i = 0;i < count;i++)
					::new (slots + i) slot();
			}
			~sscope() {
				for (int i = 0;i < count;i++)
					slots[i].~slot();
			}
		};

		template<class MF>
		OP add_arith(parsectx&ctx, OP lop, OP rop, const MF & mf) {
			return ctx.add_op([lop = std::move(lop), rop = std::move(rop), mf](slot& rv, scope *scope) {
				slot lval;
				lop->invoke(lval, scope);
				slot rval;
				rop->invoke(rval, scope);
				double dlv = lval.val()->to_num(), drv = rval.val()->to_num();
				rv.set<double>(mf(dlv, drv));
				return STATUS::OK;
			});
		}
		template<class RC, class EC>
		OP add_cmp(parsectx&ctx, OP lop, OP rop, const RC & rcmp, const EC & ecmp) {
			return ctx.add_op([lop = std::move(lop), rop = std::move(rop), rcmp, ecmp](slot& rv, scope *scope) {
				slot lval;
				lop->invoke(lval, scope);
				slot rval;
				rop->invoke(rval, scope);

				std::string * ls = lval.val()->get<std::string>(), *rs = rval.val()->get<std::string>();
				bool * lb = lval.val()->get<bool>(), *rb = rval.val()->get<bool>();
				if (ls && rs) {
					rv.set<bool>(rcmp(*ls, *rs));
				} else if (lb && rb) {
					// do a truthiness comparasion
					rv.set<bool>(ecmp(lval.val()->truthy(), rval.val()->truthy()));
				} else {
					double lnum = lval.val()->to_num(), rnum = rval.val()->to_num();
					rv.set<bool>(rcmp(lnum, rnum));
				}
				return STATUS::OK;
			});
		}

		std::shared_ptr<nfun> parse_fun(parsectx & ctx, bool statement) {

			auto argscope = ctx.enter_scope(); // this enters the arg-scope.

			int argcount = 1;
			argscope->names.emplace_back(T_THIS, "this"); // dummy first arg

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
			return std::make_shared<nfun>(id, argcount, argscope, std::move(body));
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
					val = ctx.add_op([cont = std::move(cont), val = std::move(val)](slot& rest, scope *scope) {
						// TODO: deref from cont/val
						abort();
						return STATUS::OK;
					});
					cont = nullptr;
				}
			};
			if (tt >= USER) {
				auto info = ctx.access(tt);
				if (info.first == -1) {
					auto slotptr = &(global[tt]);
					val = ctx.add_op([this, tt, slotptr](slot& res, scope *scope) {
						(*slotptr)->copyto(res);
						return STATUS::OK;
					}, tt);
				} else {
					val = ctx.add_op([info](slot&res, scope*scope) {
						int up = info.first;
						while (up--) {
							scope = scope->parent;
						}
						scope->slots[info.second].val()->copyto(res);
						return STATUS::OK;
					}, tt);
				}
			} else if (STR == tt) {
				val = ctx.add_op([tok](slot& res, scope *scope) {
					res.set<std::string>(tok);
					return STATUS::OK;
				});
			} else if (tt == DNUM) {
				val = ctx.add_op([dnum](slot&res, scope *scope) {
					res.set<double>(dnum);
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
				if (tt < prec)
					break;
				if (ADD == tt) {
					lexeat(ctx);
					deref(); // deref any references since we only want a value
					OP sub = expr(tt, ctx, MUL); // TODO: change precence req
					val = ctx.add_op([val = std::move(val), sub = std::move(sub)](slot& rv, scope *scope) {
						slot lr;
						val->invoke(lr, scope);
						slot rival;
						sub->invoke(rival, scope);
						std::string * ls = lr.val()->get<std::string>(), *rs = rival.val()->get<std::string>();
						if (ls || rs) {
							std::string sl = lr.val()->to_string(), sv = rival.val()->to_string();
							std::string out = sl + sv;
							rv.set<std::string>(out);
						} else {
							double dlv = lr.val()->to_num(), drv = rival.val()->to_num();
							rv.set<double>(dlv + drv);
						}
						return STATUS::OK;
					});
				} else if (T_DIV == tt || T_SUB == tt || T_MOD == tt || MUL == tt) {
					int op = tt;
					lexeat(ctx);
					deref();
					OP right = expr(tt, ctx, T_SUB == op ? MUL : T_MOD + 1);
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
					OP right = expr(tt, ctx, (EQ == op || NEQ == op) ? LT : GT + 1);
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
					val = ctx.add_op([args = std::move(args)](slot& rv, scope *scope) {
						slot f;
						int sz = args.size() - 1;
						slot* slots = (slot*)alloca(sizeof(slot)*sz);
						sscope dtor(slots, sz, nullptr); // responsible for construction/destruction of the alloca'd values.
						if (args[0]) {
							abort(); // TODO: "this-invoke"
						} else {
							slots[0].set<void*>(nullptr);
							args[1]->invoke(f, scope); // fetch fun.
						}
						for (int i = 1;i < sz;i++) {
							args[i + 1]->invoke(slots[i], scope);
						}

						f.val()->invoke(rv, sz, slots);
						// destroy funvalue and args
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
					auto slotptr = &(global[fn->id.first]);
					seq.insert(seq.begin(), ctx.add_op([this, fn, slotptr](slot&ref, scope *scope) {
						(*slotptr) = std::make_unique<value<decltype(fn)>>(fn);
						ref.set<decltype(fn)>(fn);
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
				seq.push_back(ctx.add_op([rexp = std::move(rexp)](slot&rv, scope * scope){
					if (rexp) {
						rexp->invoke(rv, scope);
					} else {
						rv.set<bool>(false);
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
				seq.push_back(ctx.add_op([cond = std::move(cond), tcode = std::move(tcode), fcode = std::move(fcode)](slot&res, scope *scope) {
					slot cv;
					cond->invoke(cv, scope);
					if (cv.val()->truthy()) {
						return runblock(res, scope, tcode);
					} else {
						return runblock(res, scope, fcode);
					}
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
					global[where->first] = std::move(std::make_unique<value<T>>(std::move(o)));
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

			// with the code generated
			slot rv;
			runblock(rv, nullptr, top);
		}

		runtime() {
			toksort();
			// end of built in tokens, register ID of first user-tok
			USER = nextTok;
		}
	};

	runtime::slot::~slot() {
		((valuebase*)data)->~valuebase();
	}
	runtime::slot::slot() {
		new(data) value<bool>(false);
	}
	void runtime::nfun::invoke(slot& out, int argc, slot *args)
	{
		// create local scope if the functions scope is local.
		bool is_local = info->local;
		slot* localdata = is_local ? (slot*)alloca(sizeof(slot)*info->names.size()) : nullptr;
		sscope local(localdata, localdata?info->names.size():0, info.get());
		scope * fnscope = &local;
		if (!is_local) {
			abort(); // TODO!
		}
		// TODO: varargs
		int copyargs = std::min(argc, argcount);
		for (int i = 0;i < copyargs;i++) {
			args[i].val()->moveto(fnscope->slots[i]);
		}
		runblock(out, fnscope, top);
	}
	template<class T, class ... ARG>
	inline void runtime::slot::set(ARG&&... args) {
		((valuebase*)data)->~valuebase();
		new(data) value<T>(std::forward<ARG>(args)...);
	}
};


#endif // !INCLUDED_NANOES_HPP
