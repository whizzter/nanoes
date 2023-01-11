#include <fstream>
#include <string>
#include <sstream>
#include <iostream>

#include <direct.h>

#include "nanoes.hpp"

#include <utility>
#include <vector>
#include <functional>
#include <tuple>
#include <regex>

void do_args(int argc,char** argv,std::vector < std::tuple<std::string, bool, std::function<void(std::string)>, std::string>> & params, const std::function<void(std::string arg)> & def) {
	for (int i = 1;i < argc;) {
		std::string arg(argv[i]);
		auto cmd = std::find_if(params.begin(), params.end(), [&](auto& param) {
			return 0 == arg.find_first_of(std::get<0>(param));
		});
		if (cmd != params.end()) {
			std::string val;
			if (arg.size() == std::get<0>(*cmd).size()) {
				if (std::get<1>(*cmd) && i + 1<argc) {
					val = argv[i + 1];
					i += 2;
				} else {
					i++;
				}
			} else {
				val = arg.substr(std::get<0>(*cmd).size());
				i++;
			}
			//std::cout << "Found cmd:" << std::get<0>(*cmd) << " with param " + val << "\n";
		} else {
			//std::cout << "Found reg-val:" << arg << "\n";
			i++;
		}
	}
}

std::string readfile(std::string name) {
	std::ifstream srcfile(name);
	if (!srcfile.is_open())
		return "";
	std::stringstream sstrm;
	sstrm << srcfile.rdbuf();
	return std::move(sstrm.str());
}



namespace nanoes {
	class precompiler : public runtime {
		// interpreter parsed per-operation code-snippets.
		std::map<std::string, std::vector<std::string>> ops;

		// lines of our output code.
		std::vector<std::string> lines;

		friend class runtime;
		struct prefun;
		// compiler context hooking into the regular parser/resolver.
		struct pgctx {
			size_t sp;
			std::vector<int> globals;
			//std::shared_ptr<std::map<int, std::string>> globSymMap;
			//std::vector<std::string> fnGlobs;
			precompiler *rt;

			int opc = 0;
			std::vector<std::string> lines;
			std::vector<std::shared_ptr<prefun>> funs;
			std::vector<std::shared_ptr<std::string>> slits;

			pgctx(precompiler *irt) : rt(irt) {}

			void reset() {
				opc = 0;
				lines.clear();
				globals.clear();
				funs.clear();
				slits.clear();
			}
			int label() {
				lines.push_back("op_" + std::to_string(opc)+":");
				return opc;
			}
			// this function handles the "generic" ops that we parsed from the header.
			void addTextOp(std::string op, int arg0, int arg1, bool reladdr, double *pnum) {
				// lookup the operation
				auto& oplines = rt->ops[op];
				// TODO: replace with filter that fixes ARG0/ARG1,etc
				lines.push_back("{ const int ARG0="+ std::to_string(arg0) +"; const int ARG1="+ std::to_string(arg1) +"; // begins:"+op );
				for (std::string opl : oplines) {
					//size_t pos;
					//while (std::string::npos != (pos = opl.find("ARG0"))) {
					//	opl.replace(pos, 4, std::to_string(arg0));
					//}
					//while (std::string::npos != (pos = opl.find("ARG1"))) {
					//	opl.replace(pos, 4, std::to_string(arg1));
					//}
					lines.push_back(opl);
				}
				lines.push_back("} // ends:"+op);
			}
			void add(runtime::OPC a, int32_t b = 0, int32_t c = 0, bool reladdr = false, double *dvp = nullptr) {
				opc++;
				if (a == OPC::URETURN) {
					lines.push_back("{ nanoes::INTERNALVALUE v; v.dval=0.0; return v; }");
					return;
				} else if (a == OPC::RETURN) {
					lines.push_back("return stack[" + std::to_string(b) + "];");
					return;
				} else if (a == OPC::FGOTO) {
					lines.push_back(" if (!rt->truthy(stack[" + std::to_string(b) + "])) goto op_" + std::to_string(c) + ";");
					return;
				} else if (a == OPC::GOTO) {
					lines.push_back("goto op_" + std::to_string(c) + ";");
					return;
				} else if (a == OPC::LOADDOUBLE) {
					lines.push_back("stack[" + std::to_string(b) + "].dval=" + std::to_string(*dvp) + ";" );
					return;
				} else if (a == OPC::LOADSSYM) {
					int up = c >> 10;
					int idx = c & 0x3ff;
					std::string upIn = "";
					while (up--) {
						upIn += "->parent";
					}
					lines.push_back("stack[" + std::to_string(b) + "]=cur" + upIn + "->nslots["+ std::to_string(idx)+"];");
					return;
				}
#define MATCHOP(x) if (a==OPC::x) { addTextOp(#x,b,c,reladdr,dvp); return; }
//				MATCHOP(LOADDOUBLE);
				MATCHOP(LOADGLOBAL);
				MATCHOP(LOADNULL);
				MATCHOP(LOADLIT);
				MATCHOP(LOADSSYM);
				MATCHOP(INVOKE);
				MATCHOP(ADD);
				MATCHOP(SUB);
				MATCHOP(LT);
				throw std::runtime_error("Unhandled op " + std::to_string((int)a));
			}
			int32_t addlit(std::string str) {
				int32_t idx = slits.size();
				slits.push_back(std::make_shared<std::string>(str));
				return idx;
			}
			void addfn(std::shared_ptr<prefun> && pf) {
				funs.push_back(std::move(pf));
			}
		};

		struct prefun {
//			std::shared_ptr<std::map<int, std::string>> globSymMap;
			std::vector<std::string> fnGlobs;
			std::vector<std::shared_ptr<std::string>> slits;

			std::pair<int, std::string> id;
			int argc;

			std::string ftid;
			std::string fnid;

			std::vector<std::string> lines;
			std::vector<std::shared_ptr<prefun>> fns;
			std::shared_ptr<runtime::scopeinfo> info;
			prefun(precompiler *irt, const std::pair<int, std::string> &in_id, int iargc, std::shared_ptr<runtime::scopeinfo> in_info, pgctx &ictx, std::string && src ) : id(in_id),lines(std::move(ictx.lines)) , fns(std::move(ictx.funs)), argc(iargc),info(std::move(in_info))
			{
				for (auto gi : ictx.globals) {
					auto globName = std::find_if(irt->toks.begin(), irt->toks.end(), [gi](auto& item) { return item.first == gi; })->second;
					fnGlobs.push_back(globName);
				}
				this->slits = ictx.slits;
				ftid = "FTPL" + std::to_string(uint32_t(id.first));
				fnid = "FN" + std::to_string(uint32_t(id.first));
			}
		};

		std::vector<std::shared_ptr<prefun>> topfuns;

		void dumpfn(std::ofstream& of, std::shared_ptr<prefun>&fn) {
			for (auto &sf : fn->fns) {
				dumpfn(of,sf);
			}
			of << "struct "<< fn->ftid << " { \n";
			of << "\truntime *rt;\n";
			of << "\tstd::vector<INTERNALVALUE*> globals;\n";
			of << "\tstd::vector<nesvalue> literals;\n";
			of << "\truntime::scopeinfo info;\n";
			for (auto& sf : fn->fns) {
				of << "\tstd::shared_ptr<" << sf->ftid << "> _" << sf->ftid << ";\n";
			}
			of << "\t" << fn->ftid << "(nanoes::runtime *in_rt) : rt(in_rt)\n";
			of << "\t{\n";
			of << "\t // ftpl ctor\n";
			of << "\t info.top=" << fn->info->top << ";\n";
			of << "\t //info.parent=" << fn->info->top << ";\n";
			of << "\t //info.names=" << fn->info->top << ";\n";
			of << "\t info.maxnames=" << fn->info->maxnames<< ";\n";
			of << "\t info.maxstack=" << fn->info->maxstack<< ";\n";
			of << "\t info.local=" << fn->info->local<< ";\n";
			for (auto& fg : fn->fnGlobs) {
				of << "\tglobals.push_back(&rt->global[rt->atok(\""<<fg<<"\")]);\n";
			}
			int lidx = 0;
			for (auto& slit : fn->slits) {
				of << "\t { literals.emplace_back(rt->uroot); literals["<<lidx<<"].value=rt->box(std::string(\"" << (*slit) << "\")); } \n";
				lidx++;
			}
			for (auto& sf : fn->fns) {
				//sf->info->top
				of << "\t_" << sf->ftid << " = std::make_shared<" << sf->ftid << ">(in_rt);\n";
			}
			of << "\t} // end of ftpl ctor\n";
			of << "}; // end of ftpl class\n";
			of << "struct " << fn->fnid << " : public runtime::valuebase { // "<< "dummy.." <<"\n";
			of << "std::shared_ptr<" << fn->ftid << "> bftpl;\n";
			of << "\t" << fn->fnid << "(runtime::scope *in_parent,const std::shared_ptr<" << fn->ftid << ">& in_code) : bftpl(in_code) {\n";
			of << "\t} // end of fn ctor\n";
			//of << "inline INTERNALVALUE invoke_impl() {\n";
			of << "virtual INTERNALVALUE invoke(runtime &in_rt, int argc, INTERNALVALUE *args) {\n";
			of << "\t" << fn->ftid << " *ftpl=bftpl.get();\n";
			of << "\truntime *rt=&in_rt;\n";
			if (fn->info->top) {
				of << "\truntime::scope *cur=nullptr;\n";
			} else {
				if (!fn->info->local) {
					abort();
				}
				of << "\tINTERNALVALUE localdata["<<fn->info->maxnames<<"];\n";
				of << "\truntime::scope local(localdata,&ftpl->info);\n";
				of << "\truntime::scope *cur=&local;\n";
			}
			for (auto& sf : fn->fns) {
				if (fn->info->top) {
					auto nameIdx = std::distance(fn->fnGlobs.begin(), std::find(fn->fnGlobs.begin(), fn->fnGlobs.end(), sf->id.second));
					of << "\t{ (*(ftpl->globals["<<nameIdx<<"]))=rt->box(" << sf->fnid << "(cur,ftpl->_" << sf->ftid << ")); }\n";
				} else {
					std::cout << "Found non-global?" << sf->id.second << "\n";
					abort();
				}
			}
			of << "\tINTERNALVALUE stack[" + std::to_string(fn->info->maxstack) + "];\n";
			of << "\truntime::frame cframe(rt, &cur," + std::to_string(fn->info->maxstack) + ", stack);\n";
			if (fn->argc) {
				of << "\tswitch(argc){ default:\n";
				for (int i = fn->argc - 1; i >= 0; i--) {
					of << "case " << (i + 1) << ": cur->nslots[" << i << "]=args[" << i << "];\t\n";
				}
				of << "\t}\n";
			}
			of << "\t\n";
			for (auto &l : fn->lines) {
				of << l << "\n";
			}
			of << "} // end of invoke impl \n";
			of << "virtual size_t size() { return sizeof(*this); }\n";
			of << "virtual size_t align() { return alignof(decltype(*this)); }\n";
			of << "virtual void moveto(void* dest) { ::new(dest) " << fn->fnid << "(std::move(*this)); } \n";
			of << "virtual intptr_t ty() { return id<"<<fn->fnid<<">(); }\n";
			of << " virtual void touch(runtime::toucher & to) {} // TODO\n";
			of << "virtual void * get() { return nullptr; } // no embedded objects here\n";
			of << "virtual std::string to_string() { return \" < FUNCTION... > \"; }\n";
			of << "virtual double to_num() { return 0; }\n";
			of << "virtual bool truthy() {return true;}\n";
			of << "\n";
			of << "\n";

			of << "}; // end of fn class\n";

		}
	public:
		void dump(std::string out) {
			std::ofstream of(out);

			of << "#include \"nanoes.hpp\"\n";
			of << "namespace nanoes { \n";
			of << "class precompiled { \n";
			of << "using valuebase = runtime::valuebase;\n";

			for (auto & f : topfuns) {
				dumpfn(of,f);
			}

			of << "\n";
			of << "public:\n";
			of << "static void eval(runtime &rt){\n";
			for (auto& f : topfuns) {
				auto top_name = "top_" + f->ftid;
				auto tf_name = "top_" + f->fnid;
				auto fi_name = "inst_" + f->fnid;
				of << "\tauto "<<top_name<<"=std::make_shared<"<< f->ftid <<">(&rt);\n";
				of << "\t" << f->fnid << " " << fi_name << "(nullptr," << top_name << ");\n";
				of << "\t" << fi_name << ".invoke(rt,0,nullptr);\n";
			}
			of << "\t\n";
			of << "}\n";
			of << "}; // end of precompiled\n";
			of << "}; // end of nanoes namespace\n";
		}
		bool prepare(std::string neshead) {
			std::ifstream ifile(neshead);
			if (!ifile.is_open())
				return false;

			std::string line;
			std::regex opfind("case\\s+OPC::([a-zA-Z0-9]+)");
			std::string curop;
			std::vector<std::string> *cur = nullptr;
			while (getline(ifile, line)) {
				//std::cout << "Line:[" << line << "]\n";
				std::smatch sm;
				if (std::regex_search(line, sm, opfind)) {
					//std::cout << "Found op " << sm[1] << " in line " << line << "\n";
					curop = sm[1];
					cur = &ops[sm[1]];
				} else if (cur) {
					auto ffr = line.find("continue");
					if (ffr != line.npos) {
						//std::cout << "?Ending (" << ffr << ")" << curop << " due to line :" << line << "\n";
						cur = nullptr;
						continue;
					}
					//std::cout << "?Adding to " << curop << ":" << line << "\n";
					cur->push_back(line);
				}
			}

			return true;
		}
		void compile(std::string &src,std::string & name) {
			runtime::parsectx<precompiler,pgctx, prefun> ctx(this);
			pgctx pgctx(this);
			ctx.gctx = &pgctx;
			auto fn=ctx.parse_top(src, name);

			topfuns.push_back(std::move(fn));

			//for (auto& l : fn->lines) {
			//	std::cout << "Line:[" << l << "]\n";
			//}
		}
	};
}

int main(int argc,char **argv) {
	std::string nespath = ".";

	std::vector < std::tuple<std::string,bool,std::function<void(std::string)>, std::string>> params;
	params.emplace_back("-nes", true, [&](auto arg) { nespath = arg; }, "base path to find nano es");

	do_args(argc, argv, params, [](auto def) {});

	nanoes::precompiler precomp;

	if (!precomp.prepare(nespath + "/nanoes.hpp")) {
		std::cerr << "Error, could not load :" << nespath << "/nanoes.hpp , try locating the nanoes.hpp header and give the path to it via the -nes parameter\n";
		return -1;
	}
	std::string fname = "nfib.js";
	//std::string fname = "helloworld.js";
	std::string src = readfile(fname);
	precomp.compile(src, fname);
	precomp.dump("prec.hpp");
//	char out[301] = { 0 };
//	getcwd(out, 300);
//	std::cout << out<<std::endl;
//	std::cout << src;
	return 0;
//	while (!src.eof()) {
//		std::
//		std::getline(src);
//	}
}
