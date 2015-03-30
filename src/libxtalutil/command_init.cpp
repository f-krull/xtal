#include "command.h"
#include "log.h"

#include <memory>

using namespace std;

/*----------------------------------------------------------------------------*/

namespace com {

/*----------------------------------------------------------------------------*/

//void load_cmd() {
//   if (cmd.getNumArgs() == 0) {
//      /* no argument given */
//      return;
//   }
//   source.load(cmd.getArgStr(0));
//}

/*----------------------------------------------------------------------------*/

//void zoom_cmd() {
//   const vector<Atom*> *atoms;
//
//   if (cmd.getNumArgs() == 0) {
//      atoms = source.getSelection();
//   } else {
//      atoms = source.select(cmd.getArgStr(0));
//   }
//   if (atoms->size() != 0) {
//      viewer.setCentre(Atom::getCaCentre(atoms));
//   }
//}

/*----------------------------------------------------------------------------*/
//
//void select_cmd() {
//   if (cmd.getNumArgs() == 1) {
//      source.select(cmd.getArgStr(0));
//   }
//}

/*----------------------------------------------------------------------------*/

//void show_cmd() {
//   string modelname;
//   const vector<Atom*> *atoms;
//   Model *model;
//
//   if (cmd.getNumArgs() == 0) {
//      return;
//   }
//   modelname = cmd.getArgStr(0);
//   if (cmd.getNumArgs() == 1) {
//      atoms = source.getSelection();
//   } else {
//      atoms = source.select(cmd.getArgStr(1));
//   }
//   if (modelname == "atoms") {
//      model = new AtomModel(atoms);
//   } else if (modelname == "sticks") {
//      vector<Residue*> resis;
//      Residue::getResidues(atoms, &resis);
//      Residue::calcBonds(resis);
//      model = new SticksModel(&resis);
//   } else if (modelname == "lines") {
//      vector<Residue*> resis;
//      Residue::getResidues(atoms, &resis);
//      Residue::calcBonds(resis);
//      model = new BondsModel(&resis);
//   } else {
//      return;
//   }
//   viewer.addModel(model);
//}

/*----------------------------------------------------------------------------*/

//void reset_cmd() {
//   viewer.clearModels();
//   source.clear();
//}

/*----------------------------------------------------------------------------*/

void test_cmd() {
   if (cmd.getNumArgs() < 1) {
      return;
   }
   Log::inf("test");
}

/*----------------------------------------------------------------------------*/

void run_cmd() {
   if (cmd.getNumArgs() < 1) {
      return;
   }
   cmd.readScriptfile(cmd.getArgStr(0));
}

} // end of namespace com

/*----------------------------------------------------------------------------*/

void Command::init() {
   /* dummies for internal functions */
   addCmd(NULL, "set");
   addCmd(NULL, "update");
   addCmd(com::run_cmd, "run");
   /* predefine essentials */
   addVar("logger.filename", "/dev/null");
   addVar("logger.file", "true");
   //addVar("logger.append", "no");
   addVar("logger.level", "4");
//   addVar("viewer", "no");
//   addVar("viewer.screenshotfile", "./xtal-screenshot.tga");
   /* viewer */
//   addCmd(viewer.setQuit, "viewer.quit");
//   addCmd(viewer.toggleFog, "viewer.togglefog");
//   addCmd(viewer.toggleLog, "viewer.togglelog");
//   addCmd(com::load_cmd, "load");
//   addCmd(com::zoom_cmd, "zoom");
//   addCmd(com::select_cmd, "select");
//   addCmd(com::show_cmd, "show");
   /* system */
//   addCmd(com::reset_cmd, "reset");
   /* other */
//   addCmd(tmp::test, "tmp.test");
//   addCmd(com::test_cmd, "test.orient");
//   addCmd(complexcheck::ComplexChecker::start, "compcheck");
//   addCmd(apfinder::start, "apfinder");
//   addCmd(complexclusterer::start, "compcluster");
}
