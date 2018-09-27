/****************************************************************************
**
** Copyright (C) 2009 Nokia Corporation and/or its subsidiary(-ies).
** Contact: Nokia Corporation (qt-info@nokia.com)
**
** This file is part of the examples of the Qt Toolkit.
**
** $QT_BEGIN_LICENSE:LGPL$
** Commercial Usage
** Licensees holding valid Qt Commercial licenses may use this file in
** accordance with the Qt Commercial License Agreement provided with the
** Software or, alternatively, in accordance with the terms contained in
** a written agreement between you and Nokia.
**
** GNU Lesser General Public License Usage
** Alternatively, this file may be used under the terms of the GNU Lesser
** General Public License version 2.1 as published by the Free Software
** Foundation and appearing in the file LICENSE.LGPL included in the
** packaging of this file.  Please review the following information to
** ensure the GNU Lesser General Public License version 2.1 requirements
** will be met: http://www.gnu.org/licenses/old-licenses/lgpl-2.1.html.
**
** In addition, as a special exception, Nokia gives you certain
** additional rights. These rights are described in the Nokia Qt LGPL
** Exception version 1.0, included in the file LGPL_EXCEPTION.txt in this
** package.
**
** GNU General Public License Usage
** Alternatively, this file may be used under the terms of the GNU
** General Public License version 3.0 as published by the Free Software
** Foundation and appearing in the file LICENSE.GPL included in the
** packaging of this file.  Please review the following information to
** ensure the GNU General Public License version 3.0 requirements will be
** met: http://www.gnu.org/copyleft/gpl.html.
**
** If you are unsure which license is appropriate for your use, please
** contact the sales department at http://www.qtsoftware.com/contact.
** $QT_END_LICENSE$
**
****************************************************************************/

#include <QtGui/QMessageBox>
#include <QtGui/QMessageBox>
#include <QString>
#include <QTextEdit>
#include <QFile>
#include <QFont>
#include <QMenu>
#include <QMenuBar>
#include <QFileDialog>
#include <QCoreApplication>
#include <QTextDocument>
#include <QTextCursor>
#include <QDebug>
#include <QRegExp>
#include <QDir>
#include <QFileInfo>
#include <QStringList>
 
#include "TextEdit.h"
#include "TSystem.h"

#include "StGeomHighlighter.h"

//! [0]
//__________________________________________________________________________________
TextEdit::TextEdit(QWidget *parent)
    : QMainWindow(parent)
{
    setupFileMenu();
    setupHelpMenu();
    setupEditor();

    setCentralWidget(editor);
    setWindowTitle(tr("Syntax Highlighter"));
}
//! [0]

//__________________________________________________________________________________
void TextEdit::about()
{
    QMessageBox::about(this, tr("About Syntax Highlighter"),
                tr("<p>The <b>Syntax Highlighter</b> example shows how " \
                   "to perform simple syntax highlighting by subclassing " \
                   "the QSyntaxHighlighter class and describing " \
                   "highlighting rules using regular expressions.</p>"));
}

//__________________________________________________________________________________
void TextEdit::newFile()
{
    editor->clear();
}

//__________________________________________________________________________________
void TextEdit::openFile(const QString &path)
{
    QString fileName = path;

    if (fileName.isNull()) {
       static QString filetypes = "MORTRAN Geometry  (*.g);"
                               ";STAR Geometry macro (*.C)"
                              ";";
       // create the STAR search path
       QStringList paths; paths << "." << gSystem->Getenv("STAR"); 
       QDir::setSearchPaths("geometry",paths);
       QFileInfo file("geometry:pams/geometry");
       QString defaultDir =  file.exists() ?  file.absoluteFilePath () : "";
       fileName = QFileDialog::getOpenFileName(this,
            tr("Open Geometry File"), defaultDir ,filetypes);
    }

    if (!fileName.isEmpty()) {
        QFile file(fileName);
        if (file.open(QFile::ReadOnly | QFile::Text))
            editor->setPlainText(file.readAll());
    }
}

//! [1]
//__________________________________________________________________________________
void TextEdit::setupEditor()
{
    QFont font;
    font.setFamily("Courier");
    font.setFixedPitch(true);
    font.setPointSize(12);

    editor = new QTextEdit;
    editor->setFont(font);

    highlighter = new StGeomHighlighter(editor->document());
}
//! [1]

//__________________________________________________________________________________
void TextEdit::setupFileMenu()
{
    QMenu *fileMenu = new QMenu(tr("&File"), this);
    menuBar()->addMenu(fileMenu);

    fileMenu->addAction(tr("&New"), this, SLOT(newFile()),
                        QKeySequence(tr("Ctrl+N",
                                        "File|New")));
    fileMenu->addAction(tr("&Open..."), this, SLOT(openFile()),
                        QKeySequence(tr("Ctrl+O",
                                        "File|Open")));
//    fileMenu->addAction(tr("E&xit"), qApp, SLOT(quit()),
//                        QKeySequence(tr("Ctrl+Q",
//                                        "File|Exit")));
}

//__________________________________________________________________________________
void TextEdit::setupHelpMenu()
{
    QMenu *helpMenu = new QMenu(tr("&Help"), this);
    menuBar()->addMenu(helpMenu);

    helpMenu->addAction(tr("&About"), this, SLOT(about()));
//    helpMenu->addAction(tr("About &Qt"), qApp, SLOT(aboutQt()));
}
//__________________________________________________________________________________
void TextEdit::load( const QString &f )
{
   if ( QFile::exists( f ) ) {
    QFile file(f);
    if (file.open(QFile::ReadOnly | QFile::Text)) {
        editor->setPlainText(file.readAll());
    } else {
         QMessageBox::critical(
                this,
                tr("Open failed"),
                tr("Could not open file for reading: %1").arg(QCoreApplication::translate("QFile",file.errorString()) )
                );
      }
   }
}

//__________________________________________________________________________________
void TextEdit::findBlock(const QString & expr)
{
   QRegExp exp(QString("^Block\\s+%1").arg(expr),Qt::CaseInsensitive);
   if (QTextDocument *d = editor->document() ) // currentEditor())
   {
      QTextCursor block(d);
      while (!block.isNull() ) {
         block = d->find(exp,block, QTextDocument::FindWholeWords);
         if ( ! block.isNull() ) 
         {
            QTextCursor volumeName = d->find(QRegExp("^EndBlock",Qt::CaseInsensitive),block,QTextDocument::FindWholeWords);
            if (! volumeName.isNull() )  {
               block.setPosition(volumeName.position(), QTextCursor::KeepAnchor);
               editor->setTextCursor(block);
               editor->ensureCursorVisible();
               break;
            }
         }
         break;
      }
   }
}
