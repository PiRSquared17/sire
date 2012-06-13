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

#include "Conspire/GUI/configdocument.h"
#include "Conspire/GUI/configview.h"
#include "Conspire/GUI/optionscommand.h"
#include "Conspire/GUI/optionspage.h"
#include "Conspire/GUI/exceptionpage.h"
#include "Conspire/GUI/submitpage.h"
#include "Conspire/GUI/button.h"
#include "Conspire/GUI/widgetrack.h"

#include "Conspire/exceptions.h"

#include <QUndoStack>
#include <QPainter>
#include <QGraphicsScene>
#include <QImage>

#include <QGraphicsPixmapItem>
#include <QGraphicsLinearLayout>
#include <QGraphicsBlurEffect>

using namespace Conspire;

/** Construct an empty document */
ConfigDocument::ConfigDocument(Page *parent) : Page(parent)
{
    build();
}

/** Construct a document to view and edit 'options' */
ConfigDocument::ConfigDocument(Options options, Page *parent)
               : Page(parent)
{
    build();
    this->setOptions(options);
}

/** Destructor */
ConfigDocument::~ConfigDocument()
{
    //MAYBE MAKE SURE THAT THE DOCUMENT HAS BEEN SAVED
}

/** Build the widget */
void ConfigDocument::build()
{
    disabled_view = 0;
    menu_visible = false;
    undo_stack = new QUndoStack(this);
    connect(undo_stack, SIGNAL(canUndoChanged(bool)), 
            this, SIGNAL(canUndoChanged(bool)));
    connect(undo_stack, SIGNAL(canRedoChanged(bool)), 
            this, SIGNAL(canRedoChanged(bool)));

    top_view = new PageView(false, this);
    top_view->setTitleVisible(false);
    
    view = new ConfigView(true, this);
    
    view->setTitle( "ConfigView" );
    view->setDescription( "Description of ConfigView" );
    
    connect(view, SIGNAL(add(QString)), this, SLOT(add(QString)));
    connect(view, SIGNAL(remove(QString)), this, SLOT(remove(QString)));
    connect(view, SIGNAL(update(QString,Obj)), this, SLOT(update(QString,Obj)));
    
    connect(view, SIGNAL(canBackChanged(bool)), this, SIGNAL(canBackChanged(bool)));
    connect(view, SIGNAL(canForwardChanged(bool)),
            this, SIGNAL(canForwardChanged(bool)));

    button = new Button(Conspire::tr("Submit"), this);
    connect(button, SIGNAL(clicked()), this, SLOT(submit()));

    button->setGeometry(0, this->geometry().height()-60,
                        this->geometry().width(), 60);

    top_view->setGeometry(0, 0, this->geometry().width(), this->geometry().height()-60);
    top_view->pushed(view);
}

QTransform getZoom(int w, int h, double gw, double gh, double scale)
{
    double width = (scale*gw) / double(w);
    double height = (scale*gh) / double(h);

    conspireDebug() << w << h << gw << gh << scale << width << height;

    return QTransform( width, 0.0, 0.0, height, 0.5*scale*gw, 0.5*scale*gh );
}

QTransform getZoom(const QPixmap &pixmap, const QRectF &geometry, double scale)
{
    return getZoom(pixmap.width(), pixmap.height(), 
                   geometry.width(), geometry.height(), scale);
}

/** Change the visibility of the menu */
void ConfigDocument::setMenuVisible(bool visible)
{}

void ConfigDocument::toggleMenuVisible()
{}

void ConfigDocument::resizeEvent(QGraphicsSceneResizeEvent *e)
{
    Page::resizeEvent(e);

    button->setGeometry(0, this->geometry().height()-60,
                        this->geometry().width(), 60);

    top_view->setGeometry(0, 0, this->geometry().width(), this->geometry().height()-60);
}

void ConfigDocument::moveEvent(QGraphicsSceneMoveEvent *e)
{
    Page::moveEvent(e);

    button->setGeometry(0, this->geometry().height()-60,
                        this->geometry().width(), 60);

    top_view->setGeometry(0, 0, this->geometry().width(), this->geometry().height()-60);
}

void ConfigDocument::paint(QPainter *painter, 
                           const QStyleOptionGraphicsItem *option, 
                           QWidget *widget)
{
    Page::paint(painter, option, widget);
}

/** Return the current state of the options being viewed and edited */
Options ConfigDocument::options() const
{
    return opts;
}

/** Add the new option with key 'full_key' to the document */
void ConfigDocument::add(QString full_key)
{
    Options old_state = opts;

    try
    {
        conspireDebug() << "undo_stack->push( OptionsAddCommand(" << full_key << ") )";
        undo_stack->push( new OptionsAddCommand(this,full_key) );
    }
    catch(const Conspire::Exception &e)
    {
        conspireDebug() << "error!";
        view->pushed( PagePointer( new ExceptionPage(
                Conspire::tr("Something went wrong when you tried to add "
                             "a new Option called \"%1\"").arg(full_key), e, view) ) );

        setOptions(old_state);
    }
}

/** Remove the option with key 'full_key' from the document */
void ConfigDocument::remove(QString full_key)
{
    Options old_state = opts;

    try
    {
        undo_stack->push( new OptionsRemoveCommand(this,full_key) );
    }
    catch(const Conspire::Exception &e)
    {
        view->pushed( PagePointer( new ExceptionPage(
                Conspire::tr("Something went wrong when you tried to remove "
                             "an Option called \"%1\"").arg(full_key), e, view) ) );
    
        setOptions(old_state);
    }
}

/** Update the option with key "full_key" to have the value "new_value" */
void ConfigDocument::update(QString full_key, Obj new_value)
{
    Options old_state = opts;
    
    try
    {
        undo_stack->push( new OptionsUpdateCommand(this,full_key,
                                                   new_value.asA<Value>()) );
    }
    catch(const Conspire::Exception &e)
    {
        view->pushed( PagePointer( new ExceptionPage(
                Conspire::tr("Something went wrong when you tried to set the "
                             "value of the Option called \"%1\" to \"%2\".")
                                .arg(full_key, new_value.toString()), e, view ) ) );
                                
        setOptions(old_state);
    }
}

/** Submit the current document for processing */
void ConfigDocument::submit()
{
    button->disconnect();
    button->setText(Conspire::tr("Cancel"));
    top_view->pushed( PagePointer( new SubmitPage(opts,"pmemd",view) ) );
    
    connect(button, SIGNAL(clicked()), this, SLOT(cancel()));

    QGraphicsItem::update(this->geometry());
}

/** Cancel the submission or running job */
void ConfigDocument::cancel()
{
    button->disconnect();
    button->setText(Conspire::tr("Submit"));
    top_view->popped(false);
    
    connect(button, SIGNAL(clicked()), this, SLOT(submit()));
    
    QGraphicsItem::update(this->geometry());
}

/** Undo the last operation */
void ConfigDocument::undo()
{
    undo_stack->undo();
}

/** Redo the last undone operation */
void ConfigDocument::redo()
{
    undo_stack->redo();
}

/** Move forwards to the next view in the current tab */
void ConfigDocument::forward()
{
    view->forward();
}

/** Move backwards to the last view in the current tab */
void ConfigDocument::back()
{
    view->back();
}

/** Move to the "home" view, clearing the view history if "clear_history" is true */
void ConfigDocument::home(bool clear_history)
{
    view->home(clear_history);
}

/** Set the options object that is edited by this view */
void ConfigDocument::setOptions(Options options)
{
    if (view->isEmpty())
    {
        view->pushed( PagePointer(new OptionsPage(options,view)) );
    }
    else
    {
        try
        {
            view->reread(options);
        }
        catch(const Conspire::Exception &e)
        {
            view->pushed( PagePointer(new OptionsPage(options,view)) );
        }
    }
    
    opts = options;
    
    QGraphicsWidget::update();
}
