////////////////////////////////////////////////////////////////////////
/// COPYRIGHT NOTICE
/// COPYRIGHT (c) 2018, 金小海
/// All rights reserved.
///
/// @file    HelloWorldMaker.h
/// @version 1.0
/// @author  jinxiaohai <jinxiaohaijin@outlook.com>
/// @date    Tue Jul 10 14:04:06 2018
///
/// @brief   BNL STAR code的hello world级别的代码
///
/// 修订说明:最初版本
////////////////////////////////////////////////////////////////////////
#ifndef HelloWorldMaker_def
#define HelloWorldMaker_def

#include "StMaker.h"

/// \brief HelloWorldMaker类
///
/// HelloWorldMaker类继承于虚类StMaker，该code用于入门学习。
class HelloWorldMaker : public StMaker {
 private:
  /// 事件的处理过程
  ULong_t mEventsProcessed;

 protected:

 public:
  ///
  /// \brief 构造函数
  HelloWorldMaker();

  /// \brief 析构函数
  virtual ~HelloWorldMaker();

  /// \brief 初始化分析的工具，该函数只执行一次。
  Int_t Init();

  /// \brief 在每一个事件中的分析，该make是主要的分析。
  Int_t Make();

  /// \brief 完成分析，关闭文件，执行清理工作。
  Int_t Finish();

  /// Macro for CINT compatability
  ClassDef(HelloWorldMaker, 1)
};

#endif
