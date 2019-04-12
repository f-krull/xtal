#include <sstream>
#include <fstream>
#include <cstdlib>
#include <algorithm>
#include <iostream>
#include <string.h>
#include "log.h"
#include "command.h"

using namespace std;

#define VAR_NOT_FOUND "var_not_found"

/*----------------------------------------------------------------------------*/

Command cmd;

/*----------------------------------------------------------------------------*/

Command::Command() {
   /* this is the return value if a variable was not found */
   addVar(VAR_NOT_FOUND, "DEFINED_NO_VALUE");
}

/*----------------------------------------------------------------------------*/

bool Command::existsCmd(const string &name) {
   map<string, void(*)()>::iterator it;

   it = cmds.find(name);
   return !(it == cmds.end());
}

/*----------------------------------------------------------------------------*/

bool Command::executeCmd(const string &name) {
   map<string, void(*)()>::iterator it;

   it = cmds.find(name);
   if (it == cmds.end()) {
      return false;
   }
   it->second();
   return true;
}

/*----------------------------------------------------------------------------*/

void Command::addCmd(void(*function)(), const string &name) {
   if (existsCmd(name)) {
      Log::war("function \"%s\"already exists", name.c_str());
   } else {
      list<string>::iterator it;
      it = lower_bound(cmdnames.begin(), cmdnames.end(), name);
      cmdnames.insert(it, name);
      cmds[name] = function;
   }
}

/*----------------------------------------------------------------------------*/

//string Command::parse(const string &line) {
//   string command; /* function or variable */
//   string buffer;
//
//   stringstream s(line);
//   s >> command;
//   while (s >> buffer) {
//      args.push_back(buffer);
//   }
//   return command;
//}


/*----------------------------------------------------------------------------*/

string Command::parse(const string &l) {
   string command; /* function */
   string buffer;
   int p, q;

   q = -1;
   string line = l + " ";
   while ((unsigned int) (q + 1) < line.size()) {
      /* skip spaces */
      while (((unsigned int) q + 2 < line.size()) && (line.at(q + 1) == ' ')) {
         q++;
      }
      /* do we start with a '"'? */
      string endchars = " \n";
      if (line.at(q + 1) == '\"') {
         endchars.at(0) = '\"';
         q += 1;
      } else {
         endchars.at(0) = ' ';
      }
      /* search for end of argument */
      p = q;
      q = line.find_first_of(endchars, p + 1);
      if (q == -1) {
         break;
      }
      /* store result as command or argument */
      if (command == "") {
         command = line.substr(p + 1, q - p - 1);
      } else {
         args.push_back(line.substr(p + 1, q - p - 1));
      }
   }
   return command;
}

/*----------------------------------------------------------------------------*/

void Command::execute(const string &line) {
   string cmd;

   args.clear();
   cmd = parse(line);
   if (existsCmd(cmd) == false) {
      Log::war("function \"%s\" does not exist", cmd.c_str());
      vector<string> options;
      getCmdNames(options, cmd);
      for (unsigned int i = 0; i < options.size(); i++) {
         cout << options[i] << "\n";
      }
   } else {
      /* is internal function? */
      if (cmd == "set") {
         set();
      } else if (cmd == "update") {
         update();
      } else {
         /* call function */
         executeCmd(cmd);
      }
   }
}

/*----------------------------------------------------------------------------*/

string Command::getArgStr(unsigned int i) {
   if (i >= args.size()) {
      Log::err("argument %u not found", i);
      return "";
   }
   return args.at(i);
}

/*----------------------------------------------------------------------------*/

int Command::getArgInt(unsigned int i) {
   if (i >= args.size()) {
      Log::err("argument %u not found", i);
      return 0;
   }
   return atoi(args.at(i).c_str());
}

/*----------------------------------------------------------------------------*/

float Command::getArgFloat(unsigned int i) {
   if (i >= args.size()) {
      Log::err("argument %u not found", i);
      return 0;
   }
   return atof(args.at(i).c_str());
}

/*----------------------------------------------------------------------------*/

unsigned int Command::getNumArgs() {
   return args.size();
}

/*----------------------------------------------------------------------------*/

void Command::setArgs(const vector<string> &a) {
   this->args.clear();
   this->args.assign(a.begin(), a.end());
}

/*----------------------------------------------------------------------------*/

void Command::getCmdNames(vector<string> &lst, const string &pref) {
   list<string>::iterator it;

   if (pref == "") {
      it = cmdnames.begin();
      while (it != cmdnames.end()) {
         lst.push_back(*it);
         it++;
      }
   } else {
      list<string>::iterator it, start, end;
      unsigned int l;
      start = lower_bound(cmdnames.begin(), cmdnames.end(), pref);
      it = start;
      if (it == cmdnames.end()) {
         l = 0;
      } else {
         l = min((*it).size(), pref.size());
      }
      while ((it != cmdnames.end()) && (pref.compare(0, l, (*it), 0, l)) == 0) {
         lst.push_back(*it);
         it++;
      }
   }
}

/*----------------------------------------------------------------------------*/

bool Command::existsVar(const string &name) {
   map<string, Variable>::iterator it;

   it = vars.find(name);
   return !(it == vars.end());
}

/*----------------------------------------------------------------------------*/

void Command::setVar(const string &name, const string &value) {
   map<string, Variable>::iterator it;

   it = vars.find(name);
   if (it == vars.end()) {
      addVar(name, value);
   } else {
      it->second.set(value);
   }
}

/*----------------------------------------------------------------------------*/

void Command::updateVar(const string &name, const string &value) {
   map<string, Variable>::iterator it;

   it = vars.find(name);
   if (it == vars.end()) {
      Log::err("variable \"%s\" not found. cannot update to \"%s\"",
            name.c_str(), value.c_str());
   } else {
      it->second.set(value);
   }
}

/*----------------------------------------------------------------------------*/

void Command::addVar(const string &name, const string &value) {
   if (existsVar(name)) {
      Log::dbg("variable \"%s\" already exists", name.c_str());
   } else {
      list<string>::iterator it;
      it = lower_bound(varnames.begin(), varnames.end(), name);
      varnames.insert(it, name);
      vars.insert(pair<string, Variable> (name, Variable(value)));
   }
}

/*----------------------------------------------------------------------------*/

void Command::getVarNames(vector<string> &lst, const string &prefix) {
   list<string>::iterator it;

   if (prefix == "") {
      it = varnames.begin();
      while (it != varnames.end()) {
         lst.push_back(*it);
         it++;
      }
   } else {
      list<string>::iterator it, start, end;
      unsigned int l;
      start = lower_bound(varnames.begin(), varnames.end(), prefix);
      it = start;
      if (it == varnames.end()) {
         l = 0;
      } else {
         l = min((*it).size(), prefix.size());
      }
      while ((it != varnames.end()) && (prefix.compare(0, l, (*it), 0, l)) == 0) {
         lst.push_back(*it);
         it++;
      }
   }
}

/*----------------------------------------------------------------------------*/

Variable & Command::var(const string &name) {
   map<string, Variable>::iterator it;

   it = vars.find(name);
   if (it == vars.end()) {
      Log::err("variable \"%s\" not found", name.c_str());
      return vars.find(VAR_NOT_FOUND)->second;
   }
   return it->second;
}

/*----------------------------------------------------------------------------*/

void Command::set() {
   if (this->getNumArgs() == 0) {
      /* print call variables */
      vector<string> vars;
      getVarNames(vars, "");
      for (unsigned int i = 0; i < vars.size(); i++) {
         cout << vars[i] << "\n";
      }
   } else if (this->getNumArgs() == 1) {
      /* print value of variable */
      Log::inf("value of variable \"%s\" is \"%s\"", getArgStr(0).c_str(), var(
            getArgStr(0)).str().c_str());
   } else {
      /* set value of variable */
      setVar(getArgStr(0), getArgStr(1));
   }
}

/*----------------------------------------------------------------------------*/

void Command::update() {
   if (this->getNumArgs() == 0) {
      /* print call variables */
      vector<string> vars;
      getVarNames(vars, "");
      for (unsigned int i = 0; i < vars.size(); i++) {
         cout << vars[i] << "\n";
      }
   } else if (this->getNumArgs() == 1) {
      /* print value of variable */
      Log::inf("value of variable \"%s\" is \"%s\"", getArgStr(0).c_str(), var(
            getArgStr(0)).str().c_str());
   } else {
      /* set value of variable */
      updateVar(getArgStr(0), getArgStr(1));
   }
}

/*----------------------------------------------------------------------------*/

void Command::logVars() {
   vector<string> vars;
   getVarNames(vars, "");
   for (unsigned int i = 0; i < vars.size(); i++) {
      Log::dbg("var \"%s\" \"%s\"", vars[i].c_str(),
            var(vars[i]).str().c_str());
   }
}

/*----------------------------------------------------------------------------*/

void Command::readCmdlineSettings(int argc, char **argv) {
   vector<string> args;

   /* search for settings */
   for (int32_t i = 1; i < argc; i++) {
      if (strlen(argv[i]) < 1) {
         continue;
      }
      /* is argument a config parameter? */
      if (argv[i][0] == '-' && strchr(argv[i], '=') != NULL) {
         string c = string(&argv[i][1]);
         //c.erase(c.begin());
         replace(c.begin(), c.end(), '=', ' ');
         cmd.execute("update " + c);
      } else {
         args.push_back(string(argv[i]));
      }
   }
   setArgs(args);
}

/*----------------------------------------------------------------------------*/

void Command::readScriptfile(const string &filename) {
   ifstream infile;
   string line;

   infile.open(filename.c_str(), ios::in);
   if (!infile) {
      Log::err("cannot open script file %s", filename.c_str());
      return;
   }
   Log::dbg("executing script file %s", filename.c_str());
   while (getline(infile, line)) {
      if (line.size() < 1) {
         continue;
      }
      if (line.at(0) == '#') {
         /* skip comments */
         continue;
      }
      if ((line.substr(0, 4) == "quit") || (cin.eof())) {
         break;
      }
      cmd.execute(line);
   }
   infile.close();
}
