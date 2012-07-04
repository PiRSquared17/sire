/********************************************\
  *
  *  Conspire
  *
  *  Copyright (C) 2012  Christopher Woods
  *
  *  This program is free software; you can redistribute it and/or modify
  *  it under the terms of the GNU General Public License as published by
  *  the Free Software Foundation; either version 2 of the License, or
  *  (at your option) any later version.
  *
  *  This program is distributed in the hope that it will be useful,
  *  but WITHOUT ANY WARRANTY; without even the implied warranty of
  *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  *  GNU General Public License for more details.
  *
  *  You should have received a copy of the GNU General Public License
  *  along with this program; if not, write to the Free Software
  *  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
  *
  *  For full details of the license please see the COPYING file
  *  that should have come with this distribution.
  *
  *  You can contact the authors via the developer's mailing list
  *  at http://siremol.org
  *
\*********************************************/

#include "Conspire/GUI/submitpage.h"
#include "Conspire/GUI/widgetrack.h"
#include "Conspire/GUI/button.h"
#include "Conspire/GUI/exceptionpage.h"
#include "Conspire/GUI/widgetstack.h"

#include "Conspire/option.h"
#include "Conspire/values.h"

#include "Conspire/conspire.hpp"

#include <QCoreApplication>

#include <QPainter>
#include <QRectF>

#include <QProcess>
#include <QProcessEnvironment>
#include <QTextStream>
#include <QFile>
#include <QDir>
#include <QFileInfo>
#include <QTemporaryFile>
#include <QRegExp>
#include <QComboBox>
#include <QProgressBar>
#include <QLabel>
#include <QPixmap>

#include <QGraphicsLinearLayout>
#include <QGraphicsProxyWidget>

using namespace Conspire;

static QString install_dir 
                = "/Users/chris/Work/SoftwareProject/Conspire/source/job_classes";

static QString broker = "ssi-amrmmhd.epcc.ed.ac.uk";
//static QString broker = "127.0.0.1";

void SubmitPage::build()
{
    job_id = -1;
    output_name = "results.tar.gz";
    draw_file_progress_bar = 0;
    num_bytes_expected = 0;
    num_bytes_transferred = 0;

    this->setFlag(QGraphicsItem::ItemStacksBehindParent, true);

    WidgetRack *rack = new WidgetRack(this);
    
    QGraphicsLinearLayout *l = new QGraphicsLinearLayout(::Qt::Vertical, this);
    
    this->setLayout(l);
    
    l->insertItem(0,rack);

    stack = new WidgetStack(this);
    
    WidgetRack *sub_rack = new WidgetRack(this);
    
    QLabel *label = new QLabel(Conspire::tr("When would you like the job to finish?"));
    sub_rack->addWidget(label);
    
    QComboBox *box = new QComboBox();
    box->addItem(Conspire::tr("As soon as possible!"));
    box->addItem(Conspire::tr("As soon as is practical."));
    box->addItem(Conspire::tr("Whenever is cheapest."));
    
    QGraphicsProxyWidget *box_proxy = new QGraphicsProxyWidget(this);
    box_proxy->setWidget(box);
    box_proxy->setZValue(100);
    
    sub_rack->addWidget(box_proxy);
    
    label = new QLabel(Conspire::tr("How energy efficient would you like to be?"));
    sub_rack->addWidget(label);
    
    box = new QComboBox();
    box->addItem(Conspire::tr("Use as little electricity as possible!"));
    box->addItem(Conspire::tr("Use as little electricity as is practical."));
    box->addItem(Conspire::tr("Burn as much electricity as is needed!"));

    box_proxy = new QGraphicsProxyWidget(this);
    box_proxy->setWidget(box);
    box_proxy->setZValue(100);

    sub_rack->addWidget(box_proxy);
    sub_rack->setSizePolicy( QSizePolicy(QSizePolicy::Expanding,QSizePolicy::Expanding) );
    sub_rack->setZValue(100);
    
    stack->addWidget(sub_rack);
    
    QLabel *submitting_label = new QLabel("SUBMITTING THE JOB...");
    
    submitting_label->setScaledContents(true);
    submitting_label->setPixmap( QPixmap(":data.jpg") );
    submitting_label->setSizePolicy( QSizePolicy(QSizePolicy::Expanding,
                                                 QSizePolicy::Expanding) );
    submitting_label->setMinimumSize( QSize(10,10) );
    
    stack->addWidget(submitting_label);
    
    QLabel *cloud_label = new QLabel("CLOUD");
    
    cloud_label->setScaledContents(true);
    cloud_label->setPixmap( QPixmap(":cloud.jpg") );
    cloud_label->setSizePolicy( QSizePolicy(QSizePolicy::Expanding,QSizePolicy::Expanding) );
    cloud_label->setMinimumSize( QSize(10,10) );
    
    stack->addWidget(cloud_label);
    
    QLabel *emerald_label = new QLabel("EMERALD");
    
    emerald_label->setScaledContents(true);
    emerald_label->setPixmap( QPixmap(":cluster.jpg") );
    emerald_label->setSizePolicy( QSizePolicy(QSizePolicy::Expanding,QSizePolicy::Expanding) );
    emerald_label->setMinimumSize( QSize(10,10) );
    
    stack->addWidget(emerald_label);
    
    rack->addWidget(stack);
    
    status_label = new QLabel(Conspire::tr("Click \"Start Job\" to run the job..."));
    rack->addWidget(status_label);

    progress_bar = new QProgressBar();
    progress_bar->setMinimum(0);
    progress_bar->setMaximum(100);
    rack->addWidget(progress_bar);
    
    button = new Button(Conspire::tr("Start Job"), this);
    rack->addWidget(button);
    connect(button, SIGNAL(clicked()), this, SLOT(submit()));
}

/** Constructor */
SubmitPage::SubmitPage(Page *parent) : Page(parent)
{
    build();
}

/** Construct, passing in the options used to submit the job */
SubmitPage::SubmitPage(Options options, QString clas, Page *parent) 
           : Page(parent), opts(options), job_class(clas)
{
    build();
}

/** Destructor */
SubmitPage::~SubmitPage()
{}

void SubmitPage::resizeEvent(QGraphicsSceneResizeEvent *e)
{
    Page::resizeEvent(e);
}

void SubmitPage::moveEvent(QGraphicsSceneMoveEvent *e)
{
    Page::moveEvent(e);
}

void SubmitPage::paint(QPainter *painter, 
                       const QStyleOptionGraphicsItem *option, 
                       QWidget *widget)
{
    Page::paint(painter, option, widget);
    
    if (draw_file_progress_bar)
    {
        painter->setRenderHint(QPainter::Antialiasing);
        painter->setRenderHint(QPainter::TextAntialiasing);

        painter->setOpacity(0.5);
        
        QRectF rect( 10, 10, this->geometry().width()-20, 
                     this->geometry().height()-20 );
                     
        double prog = (1.0*num_bytes_transferred) / num_bytes_expected;
        
        int span = int(360*16*prog);
        
        painter->drawPie(rect, -int(0.5*span), span);
    }
}

void SubmitPage::allUpdate()
{
    status_label->update();
    progress_bar->update();
    QCoreApplication::processEvents();
}

void SubmitPage::submit()
{
    if (job_class.isEmpty())
        return;

    try
    {

    button->disconnect();
    button->hide();
    
    status_label->setText(Conspire::tr("Creating a temporary directory..."));
    progress_bar->setValue(5);

    allUpdate();

    conspireDebug() << "SUBMIT THE CALCULATION";

    stack->switchTo(1);

    //make a temporary directory in which to stage the script
    QString tmpdir;
    {
        //use QTemporaryFile to get a unique filename
        QTemporaryFile f;
        
        if (f.open())
            tmpdir = QFileInfo(f.fileName()).absoluteFilePath();
    }
    
    if (tmpdir.isEmpty())
        throw Siren::unavailable_resource( 
                Conspire::tr("Cannot create a temporary directory!"), CODELOC );
                
    conspireDebug() << tmpdir;
                
    QDir dir(tmpdir);
    if (not dir.mkpath("."))
        throw Siren::unavailable_resource(
                Conspire::tr("Cannot create a temporary directory! (%1)")
                        .arg(tmpdir), CODELOC );

    status_label->setText(Conspire::tr("Assembling the input files..."));
    progress_bar->setValue(10);
    allUpdate();
    
    Options submit_opts = opts;

    try
    {
        //go through the options and get all of the files - copy them into
        //the temporary directory
        QStringList all_keys = opts.keysAndIndiciesWithValue(true);
        
        double progress_value = 10;
        double progress_delta = all_keys.count() / 40.0;
        
        foreach (QString key, all_keys)
        {
            Option opt = opts[key];
            
            progress_value += progress_delta;
            progress_bar->setValue(progress_value);
            allUpdate();
            
            if (opt.value().isA<FileValue>())
            {
                conspireDebug() << opt.key() << "is a file" << opt.value().toString();
                
                //copy this file to the tmp directory
                QFile f(opt.value().toString());
                
                status_label->setText(Conspire::tr("Copying file \"%1\"...")
                                            .arg(opt.value().toString()));
                allUpdate();
                
                if (not f.open(QIODevice::ReadOnly))
                    emit( push( PagePointer( new ExceptionPage(
                        Conspire::tr("Error in submission! Cannot find the file "
                                     "\"%2\" which is required for option \"%1\".")
                                        .arg(opt.key(), opt.value().toString()),
                            Conspire::file_error( Conspire::tr("Cannot open the "
                                "file \"%2\" which is specified for option \"%1\".")
                                    .arg(opt.key(), opt.value().toString()), CODELOC ) 
                                        ) ) ) );

                QString short_name = QFileInfo(f).fileName();
                
                if (not f.copy(QString("%1/%2").arg(tmpdir,short_name)))
                    emit( push( PagePointer( new ExceptionPage(
                        Conspire::tr("Error in submission! Cannot copy the file "
                                     "\"%2\" which is required for option \"%1\".")
                                        .arg(opt.key(), opt.value().toString()),
                            Conspire::file_error( Conspire::tr("Cannot copy the "
                                "file \"%2\" which is specified for option \"%1\" "
                                "to the temporary directory \"%3\".")
                                    .arg(opt.key(), opt.value().toString(), tmpdir), 
                                        CODELOC ) ) ) ) );
            
                //update the options object so that it uses the local path for the file
                submit_opts = submit_opts.setNestedValue(key, short_name)
                                         .asA<Options>();
            }
        }
        
        status_label->setText( Conspire::tr("Creating a config file...") );
        progress_bar->setValue(50);
        allUpdate();
         
        //now write the configuration file into this directory
        QFile conf(QString("%1/config").arg(tmpdir));
        
        if (not conf.open(QIODevice::WriteOnly))
            emit( push( PagePointer( new ExceptionPage(
                Conspire::tr("Error in submission! Could not write the config file."),
                            Conspire::file_error( Conspire::tr("Cannot open the "
                                "config file \"%1/config\".")
                                    .arg(tmpdir), CODELOC ) ) ) ) );

        //write the config file
        {
            QTextStream ts(&conf);
            ts << submit_opts.toConfig(true);
            conf.close();
        }

        status_label->setText( Conspire::tr("Copying support files...") );
        progress_bar->setValue(60);
        allUpdate();

        //now copy all support files needed for this class of job
        //into the temporary file
        QDir class_dir(QString("%1/%2").arg(install_dir,job_class));
        class_dir.setFilter(QDir::Files);
        
        if (not class_dir.exists())
        {
            emit( push( PagePointer( new ExceptionPage(
                Conspire::tr("Error in submission! Could not find the job class "
                             "directory \"%1\".").arg(class_dir.absolutePath()),
                            Conspire::file_error( Conspire::tr("Cannot open the "
                                "class file directory \"%1\".")
                                    .arg(class_dir.absolutePath()), CODELOC ) ) ) ) );
        }
        
        conspireDebug() << "Reading job class directory" << class_dir.absolutePath();
        
        foreach (QString file, class_dir.entryList())
        {
            QFile f(QString("%1/%2").arg(class_dir.absolutePath(),file));
            
            QString f_name = QFileInfo(f).fileName();
            
            conspireDebug() << "Copying file" << f_name;
            
            if (not f.open(QIODevice::ReadOnly))
                throw Conspire::file_error( Conspire::tr(
                        "Cannot open class file \"%1\".")
                            .arg(file), CODELOC );
                            
            if (not f.copy(QString("%1/%2").arg(tmpdir,f_name)))
                throw Conspire::file_error( Conspire::tr(
                        "Cannot copy class file \"%1\" to temp directory \"%2\".")
                            .arg(file, tmpdir), CODELOC );
        
        }
                  
        //now tar up this whole directory into a workpacket
        {
            status_label->setText( Conspire::tr("Constructing the workpacket...") );
            progress_bar->setValue(80);
            allUpdate();
    
            conspireDebug() << "Creating the workpacket...";
            QProcess proc;
            proc.setWorkingDirectory(tmpdir);
            
            QStringList args;

            QDir d(tmpdir);
            d.setFilter(QDir::Files);

            args << "-zLcvf" << "workpacket.tgz" << d.entryList();
            conspireDebug() << "tar" << args;
            
            proc.start("tar", args);
            
            proc.waitForFinished(-1);
            
            conspireDebug() << "STDOUT" << proc.readAllStandardOutput();
            conspireDebug() << "STDERR" << proc.readAllStandardError();
            
            conspireDebug() << "...workpacket created!";
        }
        
        status_label->setText( Conspire::tr("Sending the job to the cloud...") );
        progress_bar->setValue(85);
        allUpdate();
        
        QProcessEnvironment env;
        env.insert("PARENT_NODE", broker);
        env.insert("PARENT_NODE_PORT", "10000");
        env.insert("ISSUBMIT", "TRUE");
        env.insert("WORKNAME", "workpacket");
        env.insert("CLASSNAME", job_class);
        env.insert("FILENAME", QString("%1/workpacket.tgz").arg(tmpdir));
        
        QProcess proc;
        proc.setProcessEnvironment(env);
        
        proc.start("python", 
            QStringList(QString("%1/leafhead3.py").arg(install_dir)));
        
        progress_bar->setValue(90);
        allUpdate();
        
        /*QTextStream ts(&proc);
        QStringList lines;
        
        proc.closeWriteChannel();
        proc.waitForStarted(-1);
        
        draw_file_progress_bar = true;
        num_bytes_expected = QFileInfo(QString("%1/workpacket.tgz").arg(tmpdir)).size();
        num_bytes_transferred = 0;
        
        while (proc.state() == QProcess::Running)
        {
            if (proc.waitForReadyRead(100))
            {
                if (proc.canReadLine())
                {
                    QString line = ts.readLine();
                    conspireDebug() << "OUTPUT" << line;
                    lines.append(line);
                    
                    bool ok = false;
                    int completed = line.toInt(&ok);
                    
                    if (ok)
                    {
                        num_bytes_transferred = completed;
                        this->repaint();
                    }
                }
            }
        }
        
        draw_file_progress_bar = false;*/
        
        proc.waitForFinished(-1);
        
        QByteArray out = proc.readAllStandardOutput();
        //QString out = lines.join("\n");
        QByteArray err = proc.readAllStandardError();
        
        conspireDebug() << "OUTPUT" << out;
        conspireDebug() << "ERROR" << err;
        
        conspireDebug() << "JOB SUBMITTED";
        
        if (proc.exitCode() != 0)
        {
            status_label->setText( Conspire::tr("SUBMISSION FAILURE!!!") );
            progress_bar->setValue(0);
            allUpdate();
        
            //something went wrong
            emit( push( PagePointer( new ExceptionPage("Error in submission!", 
                                            unavailable_resource("Could not submit!", 
                                                CODELOC) ) ) ) );
        }

        stack->switchTo(2);
        status_label->setText( Conspire::tr("Job is on the cloud!") );
        progress_bar->setValue(95);
        allUpdate();

        //get the JOB_ID and pass it to the job running page
        QRegExp re("ID on BROKER: (\\d+)");

        int pos = re.indexIn(out);
        
        if (pos == -1)
            throw Conspire::unavailable_resource("ERROR IN PARSING RESULT OF SUBMISSION!",
                                                    CODELOC );

        conspireDebug() << "JOB_ID ==" << re.cap(1);
        job_id = re.cap(1).toInt();
        
        status_label->setText( Conspire::tr("Submitted the job "
               "to the cloud! JOB_ID = %1.").arg(job_id) );
        progress_bar->setValue(100);
        allUpdate();
        
        button->setText( Conspire::tr("Get Status") );
        connect(button, SIGNAL(clicked()), this, SLOT(query()));
        button->show();
    }
    catch(...)
    {
        dir.rmpath(".");
        throw;
    }
    
    dir.rmpath(".");
    }
    catch(const Siren::Exception &e)
    {
        emit( push( PagePointer( new ExceptionPage(
            Conspire::tr("Error in submission!"), e) ) ) );
    }
}

void SubmitPage::query()
{
    button->hide();
    
    status_label->setText( Conspire::tr("Querying the job...") );
    progress_bar->setValue(0);
    allUpdate();
    
    //do the query
    //PARENT_NODE=127.0.0.1 PARENT_NODE_PORT=10000 ISQUERY=TRUE WKPTID=2 python ./leafhead3.py
    QProcessEnvironment env;
    env.insert("PARENT_NODE", broker);
    env.insert("PARENT_NODE_PORT", "10000");
    env.insert("ISQUERY", "TRUE");
    env.insert("WKPTID", QString::number(job_id));
    
    QProcess proc;
    proc.setProcessEnvironment(env);
    
    proc.start("python", 
        QStringList(QString("%1/leafhead3.py").arg(install_dir)));
    
    progress_bar->setValue(50);
    allUpdate();
    
    proc.waitForFinished(-1);
    
    QByteArray out = proc.readAllStandardOutput();
    QByteArray err = proc.readAllStandardError();
    
    conspireDebug() << "OUTPUT" << out;
    conspireDebug() << "ERROR" << err;

    status_label->setText( Conspire::tr("Got a response...") );
    progress_bar->setValue(50);
    allUpdate();

    //WORKPACKET ID 0 on BROKER is *******
    // where ****** is unallocated, allocated, running, completed, or non existant
    QRegExp re("WORKPACKET ID (\\d+) on BROKER is ([\\w\\s]+)");
    
    int pos = re.indexIn(out);
    
    if (pos == -1)
    {
        status_label->setText( Conspire::tr("Strange? Couldn't get job status.") );
        progress_bar->setValue(0);
        allUpdate();
    }
    else
    {
        QString status = re.cap(2).simplified();

        if (status == "running")
        {
            stack->switchTo(3);
        }
    
        status_label->setText( Conspire::tr("The job status is \"%1\".").arg(status) );
        progress_bar->setValue(100);
        allUpdate();
    
        if (status == "completed")
        {
            stack->switchTo(1);
            button->setText( Conspire::tr("Get Results") );
            button->disconnect();
            connect(button, SIGNAL(clicked()), this, SLOT(getResults()));
        }
    }
    
    button->show();
}

void SubmitPage::getResults()
{
    //button->disconnect();
    button->hide();
    
    status_label->setText( Conspire::tr("Getting the results...") );
    progress_bar->setValue(0);
    allUpdate();

    //do the query
    //PARENT_NODE=127.0.0.1 PARENT_NODE_PORT=10000 ISQUERY=TRUE WKPTID=2 python ./leafhead3.py
    QProcessEnvironment env;
    env.insert("PARENT_NODE", broker);
    env.insert("PARENT_NODE_PORT", "10000");
    env.insert("ISDOWNLOAD", "TRUE");
    env.insert("WKPTID", QString::number(job_id));
    env.insert("OUTFILENAME", output_name);

    QProcess proc;
    proc.setProcessEnvironment(env);
    
    proc.start("python", 
        QStringList(QString("%1/leafhead3.py").arg(install_dir)));
    
    progress_bar->setValue(50);
    allUpdate();
    
    proc.waitForFinished(-1);
    
    QByteArray out = proc.readAllStandardOutput();
    QByteArray err = proc.readAllStandardError();
    
    conspireDebug() << "OUTPUT" << out;
    conspireDebug() << "ERROR" << err;

    status_label->setText( Conspire::tr("Downloaded results to %1").arg(output_name) );
    progress_bar->setValue(100);

    button->show();
    allUpdate();

    QProcess proc2;
    proc2.start("open", QStringList() << output_name);

    proc2.waitForFinished(-1);
}