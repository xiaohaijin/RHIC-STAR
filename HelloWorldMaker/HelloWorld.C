#include <TROOT.h>
#include <TStyle.h>
#include <TSystem.h>

#include <StChain.h>

#include "StRoot/HelloWorldMaker/HelloWorldMaker.h"
void HelloWorld() {
  /// Execute a macro in the interpreter(在解释器中执行loadMuDst.C)
  gROOT->Macro("loadMuDst.C");
  /// 加载动态库
  gSystem->Load("HelloWorldMaker.so");

  // List of member links in the chain
  StChain* chain = new StChain;
  HelloWorldMaker* Hello = new HelloWorldMaker();

  /// 事件数
  Int_t nEvents = 10;

  /// Loop over the links in the chain \n
  /// Main base class to control chains for the different STAR "chains"
  chain->Init();

  chain->EventLoop(1, nEvents);

  chain->Finish();

  // 清理工作
  delete chain;
}
