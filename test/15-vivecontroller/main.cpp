#include <ViveController/vivecontroller.h>

#include <Core/array.h>


const char *USAGE =
    "\nTest of low-level (without bot interface) ViveController interfacee"
    "\n";

int main(int argc,char **argv){
  rai::initCmdLine(argc, argv);

  cout <<USAGE <<endl;

  rai::Configuration C;
  rai::ViveController VC;

  for(;;){
      VC.pull(C);
    if(C.view(false)=='q') break;
  }

  LOG(0) <<" === bye bye ===\n used parameters:\n" <<rai::params() <<'\n';

  return 0;
}
