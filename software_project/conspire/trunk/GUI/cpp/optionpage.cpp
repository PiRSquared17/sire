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

#include "Conspire/GUI/optionpage.h"
#include "Conspire/GUI/exceptionpage.h"
#include "Conspire/GUI/slidelabel.h"

#include "Conspire/values.h"
#include "Conspire/exceptions.h"

#include <QGraphicsLinearLayout>
#include <QGraphicsProxyWidget>

#include <QLabel>
#include <QLineEdit>

using namespace Conspire;

/** Constructor */
OptionPage::OptionPage(Page *parent) : ConfigPage(parent)
{
    build();
}

/** Constructor, passing in the option to edit, and the root key of that option */
OptionPage::OptionPage(Option option, QString root_key, Page *parent)
           : ConfigPage(parent)
{
    build();
    setOption(option, root_key);
}
   
/** Destructor */
OptionPage::~OptionPage()
{}

/** Return the option being viewed and edited */
Option OptionPage::option() const
{
    return opt;
}

/** Return the root key of that option */
QString OptionPage::rootKey() const
{
    return root_key;
}

/** Return the short key of the option (the key without the root_key) */
QString OptionPage::shortKey() const
{
    if (opt.index() == 1)
        return opt.key();
    else
        return QString("%1[%2]").arg(opt.key()).arg(opt.index());
}

/** Return the full key of the option, including the root key */
QString OptionPage::fullKey() const
{
    if (root_key.isEmpty())
        return shortKey();
    else
        return QString("%1.%2").arg(root_key, shortKey());
}

/** Build the widget */
void OptionPage::build()
{
    setTitle("Unnamed OptionPage");
    setDescription("Unnamed OptionPage description");
    value_edit = new QLineEdit();
    value_edit->setAlignment( ::Qt::AlignCenter );
    value_edit->setFont( QFont("LucidaGrande", 48) );

    connect(value_edit, SIGNAL(editingFinished()), this, SLOT(editingFinished()));

    help_label = new SlideLabel(this);
    
    QGraphicsProxyWidget *edit_proxy = new QGraphicsProxyWidget(help_label);
    edit_proxy->setWidget(value_edit);
    
    help_label->setWidget(edit_proxy);
}

void OptionPage::resizeEvent(QGraphicsSceneResizeEvent *e)
{
    Page::resizeEvent(e);
    help_label->setGeometry(0, 0, this->geometry().width(), this->geometry().height());
}

void OptionPage::moveEvent(QGraphicsSceneMoveEvent *e)
{
    Page::moveEvent(e);
    help_label->setGeometry(0, 0, this->geometry().width(), this->geometry().height());
}

/** Function called when editing of the option is finished */
void OptionPage::editingFinished()
{
    conspireDebug() << "EDITING FINISHED" << value_edit->text() << old_text;

    if (value_edit->text() == old_text)
        return;

    old_text = value_edit->text();

    //try to update the option locally first. This catches parse errors
    //before they propogate up the document
    try
    {
        Obj new_val = StringValue(value_edit->text());

        Options opts(opt);
        opts.setNestedValue(shortKey(),new_val.asA<Value>());
        
        conspireDebug() << "emit( update(" << fullKey() << "," << value_edit->text() 
                        << " ) )";
        emit( update(fullKey(),new_val) );
        emit( pop() );
    }
    catch(const Conspire::Exception &e)
    {
        emit( push(PagePointer( new ExceptionPage( Conspire::tr(
                    "Something went wrong when updating the value of the option \"%1\""
                    "to \"%2\".")
                        .arg(fullKey(), value_edit->text()), e) ) ) );
                        
        value_edit->setText( opt.value().toString() );
        old_text = value_edit->text();
    }
}

/** Set the option to be edited */
void OptionPage::setOption(Option option, QString key)
{
    opt = option;
    root_key = key;
    
    this->setDescription(option.description());
    this->setTitle(fullKey());
    
    value_edit->setText(opt.value().toString());
    old_text = value_edit->text();
}

/** Slot called when the option has been updated */
void OptionPage::reread(Options options)
{
    try
    {
        Option opt = options[fullKey()];
    
        QString text = opt.value().toString();
    
        if (text != old_text)
        {
            value_edit->setText(text);
            old_text = text;
        }
    }
    catch(...)
    {
        conspireDebug() << "ERROR IN REREADING OPTION...";
    }
    
    ConfigPage::reread(options);
}
