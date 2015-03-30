#ifndef COMMAND_H_
#define COMMAND_H_

#include <string>
#include <vector>
#include <map>
#include <list>
#include "variable.h"

/*----------------------------------------------------------------------------*/

class Command {
public:
   Command();
   void init();

   /* commands */
   bool existsCmd(const std::string &name);
   void addCmd(void(*function)(), const std::string &name);
   bool executeCmd(const std::string &name);
   void getCmdNames(std::vector<std::string> &n, const std::string &pref = "");

   /* variables */
   bool existsVar(const std::string &name);
   void addVar(const std::string &name, const std::string &value);
   void setVar(const std::string &name, const std::string &value);
   void updateVar(const std::string &name, const std::string &value);
   Variable& var(const std::string &name);
   void getVarNames(std::vector<std::string> &n, const std::string &pref = "");

   /* execute string (calls function and sets arguments) */
   void execute(const std::string &line);

   /* arguments of functions / name and value of variable for "set" */
   std::string getArgStr(unsigned int i);
   int getArgInt(unsigned int i);
   float getArgFloat(unsigned int i);
   unsigned int getNumArgs();
   void setArgs(const std::vector<std::string> &args);

   /* write all defined variables to log */
   void logVars();
   /* set command line arguments as variables ("-logger.file=false") */
   void readCmdlineSettings(int argc, char **argv);
   /* execute each line of a file */
   void readScriptfile(const std::string &filename);

private:
   std::map<std::string, void(*)()> cmds;
   std::list<std::string> cmdnames;

   std::map<std::string, Variable> vars;
   std::list<std::string> varnames;

   std::vector<std::string> args;
   std::string parse(const std::string &line);
   void set();
   void update();
};

extern Command cmd;

#endif /* COMMAND_H_ */
