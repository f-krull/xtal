#include "exectemplate.h"
#include "../libxtalutil/log.h"
#include "../libxtalutil/common.h"
#include "../libxtalutil/command.h"
#include <typeinfo>
#include <stdlib.h>

/*----------------------------------------------------------------------------*/

void ExecTemplate::registerStuff() {
}

static bool init(ExecTemplate *b, int argc, char **argv) {
   cmd.init();
   b->registerStuff();
   cmd.addVar("sys.configfile", "./data/config.txt");
   if (common::exists(cmd.var("sys.configfile").str()) == true) {
      cmd.readScriptfile(cmd.var("sys.configfile").str());
   }
   cmd.readCmdlineSettings(argc, argv);
   Log::inf("build on %s %s for %ubit", __DATE__, __TIME__, sizeof(void*) * 8);
   cmd.logVars();
   return true;
}

int ExecTemplate::start(int argc, char **argv) {
   int ret;
   if (init(this, argc, argv) == false) {
      return EXIT_FAILURE;
   }
   ret = start();
   return ret;
}

int start(ExecTemplate *e, int argc, char **argv) {
   return e->start(argc, argv);
}

/*----------------------------------------------------------------------------*/

#include <mpi.h>
int MpiTemplate::start(int argc, char **argv) {
   MPI_Init (&argc, &argv);
   int ret = ExecTemplate::start(argc, argv);
   MPI_Finalize();
   return ret;
}

