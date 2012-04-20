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

#include "Conspire/GUI/optionspage.h"
#include "Conspire/GUI/addpage.h"
#include "Conspire/GUI/optionpage.h"
#include "Conspire/GUI/exceptionpage.h"

#include <QGraphicsLinearLayout>
#include <QGraphicsProxyWidget>

#include <QPushButton>
#include <QFrame>
#include <QScrollArea>
#include <QVBoxLayout>

#include <QSignalMapper>

using namespace Conspire;

/** Constructor */
OptionsPage::OptionsPage(Page *parent) : ConfigPage(parent)
{
    build();
}

/** Construct, passing in a top-level options object */
OptionsPage::OptionsPage(Options options, Page *parent)
            : ConfigPage(parent)
{
    build();
    setOptions(options);
}

/** Construct, passing in an Opitons object that is a sub-options object,
    with root key 'root_key' */
OptionsPage::OptionsPage(Options options, QString root_key, Page *parent)
            : ConfigPage(parent)
{
    build();
    setOptions(options, root_key);
}

/** Destructor */
OptionsPage::~OptionsPage()
{}

/** Return the options object being edited by this view */
Options OptionsPage::options() const
{
    return opts;
}

/** Return the root key of the options object being edited. This is
    null if this view is editing the top-level options object */
QString OptionsPage::rootKey() const
{
    return root_key;
}

/** Internal function used to rebuild the widget from the passed options */
void OptionsPage::pvt_reread(Options options)
{
    //create one button for each option...
    opts = options;
    
    QStringList keys = opts.keysAndIndiciesWithValue();

    if (keys.count() > buttons.count())
    {
        for (int i=buttons.count(); i<keys.count(); ++i)
        {
            QPushButton *b = new QPushButton(button_frame);
            button_frame->layout()->addWidget(b);
            buttons.append(b);

            connect(b, SIGNAL(clicked()), mapper, SLOT(map()));
        }
    }
    else if (keys.count() < buttons.count())
    {
        while (buttons.count() > keys.count())
        {
            QAbstractButton *b = buttons.takeLast();
            delete b;
        }
    }
    
    for (int i=0; i<keys.count(); ++i)
    {
        QAbstractButton *b = buttons[i];
        Option opt = opts.getNestedOption(keys[i]);
        
        if (opt.hasSubOptions())
        {
            b->setText( QString("%1 ----->").arg(keys[i]) );
        }
        else
        {
            b->setText( QString("%1 == %2").arg(keys[i], opt.value().toString()) );
        }
        
        b->show();
        
        mapper->removeMappings(b);
        mapper->setMapping(b, keys[i]);
    }
    
    add_button->setVisible( opts.canAddValues() );
}

/** This function is called when the options object is updated */
void OptionsPage::reread(Options options)
{
    if (not root_key.isNull())
    {
        try
        {
            Option opt = options.getNestedOption(root_key);
            options = opt.value().asA<Options>();
        }
        catch(...)
        {
            conspireDebug() << "CANNOT FIND" << root_key << "IN OPTIONS!";
            options = Options();
        }
    }
    
    pvt_reread(options);
}

/** Slot called when one of the options is clicked */
void OptionsPage::clicked(const QString &key)
{
    try
    {
        //create a new OptionsWidget to display the suboptions
        Option opt = opts.getNestedOption(key);
        
        conspireDebug() << opt.toString();
        
        if (opt.hasSubOptions())
        {
            QString root = key;
                
            if (not root_key.isEmpty())
                root.prepend(".").prepend(root_key);

            OptionsPage *page = new OptionsPage(opt.value().asA<Options>(), root);

            page->setTitle(root);
            page->setDescription(opt.description());

            emit( push(PagePointer(page)) );
        }
        else
        {
            OptionPage *page = new OptionPage(opt, root_key);
            
            if (not root_key.isEmpty())
            {
                page->setTitle( QString("%1.%2").arg(root_key, opt.key()) );
            }
            else
            {
                page->setTitle(opt.key());
            }
            
            page->setDescription(opt.description());
            
            emit( push(PagePointer(page)) );
        }
    }
    catch(const Conspire::Exception &e)
    {
        conspireDebug() << "EXCEPTION" << e.toString();
        
        PagePointer p(new ExceptionPage( Conspire::tr(
                "An error occurred when you clicked on the option \"%1\".")
                    .arg(key), e ));
                    
        emit( push(p) );
    }
}

/** Actually build the widget */
void OptionsPage::build()
{
    QFrame *page_frame = new QFrame();
    page_frame->setFrameShadow(QFrame::Plain);
    page_frame->setFrameShape(QFrame::NoFrame);
    page_frame->setLineWidth(0);
    page_frame->setMidLineWidth(0);
    page_frame->setContentsMargins(0,0,0,0);
    
    page_frame->setLayout( new QVBoxLayout(page_frame) );
    
    QScrollArea *scroller = new QScrollArea(page_frame);
    page_frame->layout()->addWidget(scroller);
    scroller->show();
    
    add_button = new QPushButton( Conspire::tr("Add..."), page_frame );
    add_button->setSizePolicy( QSizePolicy::MinimumExpanding,
                               QSizePolicy::Minimum );
    connect(add_button, SIGNAL(clicked()), this, SLOT(add()));

    page_frame->layout()->addWidget(add_button);
    add_button->hide();

    button_frame = new QFrame(scroller);
    button_frame->setFrameShadow(QFrame::Plain);
    button_frame->setFrameShape(QFrame::Box);
    button_frame->setLineWidth(1);
    button_frame->setMidLineWidth(0);
    button_frame->setContentsMargins(2,2,2,2);

    button_frame->setLayout( new QVBoxLayout(button_frame) );
    button_frame->show();
    scroller->setWidget(button_frame);
    scroller->setWidgetResizable(true);
    
    mapper = new QSignalMapper(this);
    connect(mapper, SIGNAL(mapped(const QString&)), this, SLOT(clicked(const QString&)));

    this->setPageWidget(page_frame);
}

/** Set the options and root_key used by this object */
void OptionsPage::setOptions(Options options, QString key)
{
    root_key = key;
    opts = Options();
    this->pvt_reread(options);
}

/** Called when it the user wants to add an option to this page */
void OptionsPage::add()
{
    if (opts.canAddValues())
    {
        AddPage *page = new AddPage(opts, root_key);
        emit( push(PagePointer(page)) );
    }
}
