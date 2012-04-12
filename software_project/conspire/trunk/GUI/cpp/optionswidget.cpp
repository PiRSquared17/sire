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

#include "Conspire/GUI/optionswidget.h"
#include "Conspire/GUI/addwidget.h"
#include "Conspire/GUI/optionseditview.h"
#include "Conspire/GUI/optionscommand.h"
#include "Conspire/GUI/mainbar.h"

#include <QVBoxLayout>

#include <QGraphicsView>
#include <QGraphicsScene>
#include <QGraphicsGridLayout>
#include <QGraphicsLinearLayout>
#include <QGraphicsWidget>
#include <QGraphicsProxyWidget>
#include <QUndoStack>
#include <QLabel>
#include <QPushButton>

#include <QPropertyAnimation>
#include <QParallelAnimationGroup>

using namespace Conspire;

/////////
///////// Implementation of NewOptionsWidget
/////////

/** Constructor */
NewOptionsWidget::NewOptionsWidget(QGraphicsItem *parent) : PageWidget(parent)
{
    build();
}

/** Construct, passing in a top-level options object */
NewOptionsWidget::NewOptionsWidget(Options options, QGraphicsItem *parent)
                 : PageWidget(parent)
{
    build();
    setOptions(options);
}

/** Construct, passing in an Opitons object that is a sub-options object,
    with root key 'root_key' */
NewOptionsWidget::NewOptionsWidget(Options options, QString root_key,
                                   QGraphicsItem *parent)
                 : PageWidget(parent)
{
    build();
    setOptions(options, root_key);
}

/** Destructor */
NewOptionsWidget::~NewOptionsWidget()
{}

/** Return the options object being edited by this view */
Options NewOptionsWidget::options() const
{
    return opts;
}

/** Return the root key of the options object being edited. This is
    null if this view is editing the top-level options object */
QString NewOptionsWidget::rootKey() const
{
    return root_key;
}

/** This function is called when the options object is updated */
void NewOptionsWidget::update(Options options)
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
    
    //create one button for each option...
    opts = options;
    
    QStringList keys = opts.keysAndIndiciesWithValue();
    
    //make sure that there are as many buttons in the layout as options
    QGraphicsLinearLayout *l = dynamic_cast<QGraphicsLinearLayout*>(this->layout());
    
    if (not l)
    {
        conspireDebug() << "INVALID LAYOUT!!!";
        throw Conspire::program_bug( Conspire::tr(
                "How do we have an invalid layout???"), CODELOC );
    }
    
    if (l->count() > keys.count())
    {
        for (int i = l->count()-1; i >= keys.count(); --i)
        {
            QGraphicsLayoutItem *item = l->itemAt(i);
            l->removeAt(i);
            buttons.removeAt(i);
            delete item;
        }
    }
    else if (l->count() < keys.count())
    {
        for (int i = l->count(); i < keys.count(); ++i)
        {
            QPushButton *b = new QPushButton();
            connect(b, SIGNAL(clicked()), mapper, SLOT(map()));
            buttons.append(b);
            QGraphicsProxyWidget *b_proxy = new QGraphicsProxyWidget(this);
            b_proxy->setWidget(b);
            l->insertItem(i, b_proxy);
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
        
        mapper->removeMappings(b);
        mapper->setMapping(b, keys[i]);
    }
}

/** Slot called when one of the options is clicked */
void NewOptionsWidget::clicked(const QString &key)
{
    try
    {
        //create a new OptionsWidget to display the suboptions
        conspireDebug() << "CLICKED" << key;
        Option opt = opts.getNestedOption(key);
        
        if (opt.hasSubOptions())
        {
            QString root = key;
                
            if (not root_key.isNull())
                root.prepend(".").prepend(root_key);

            PagePointer p(new NewOptionsWidget(opt.value().asA<Options>(), root));
            
            emit( push(p) );
        }
        else
        {
            PagePointer p(new OptionWidget(opt, root_key));
            
            emit( push(p) );
        }
    }
    catch(const Conspire::Exception &e)
    {
        conspireDebug() << "EXCEPTION" << e.toString();
        
        PagePointer p(new ExceptionWidget( Conspire::tr(
                "An error occurred when you clicked on the option \"%1\".")
                    .arg(key), e ));
                    
        emit( push(p) );
    }
}

/** Actually build the widget */
void NewOptionsWidget::build()
{
    QGraphicsLinearLayout *l = new QGraphicsLinearLayout( ::Qt::Vertical, this );
    this->setLayout(l);
    
    mapper = new QSignalMapper(this);
    connect(mapper, SIGNAL(mapped(const QString&)), this, SLOT(clicked(const QString&)));
}

/** Set the options and root_key used by this object */
void NewOptionsWidget::setOptions(Options options, QString key)
{
    root_key = key;
    opts = Options();
    this->update(options);
}

/////////
///////// Implementation of OptionsWidget
/////////

/** Constructor */
OptionsWidget::OptionsWidget(QWidget *parent) : QGraphicsView(parent)
{
    build();
}

/** Construct to edit the passed options */
OptionsWidget::OptionsWidget(Options options, QWidget *parent)
              : QGraphicsView(parent)
{
    build();
    setOptions(options);
}

/** Destructor */
OptionsWidget::~OptionsWidget()
{}

/** Return the current state of the options being edited */
Options OptionsWidget::options() const
{
    return opts;
}

/** Set the options object to be edited */
void OptionsWidget::setOptions(Options options)
{
    try
    {
        //top_view->setOptions(options);
        top_view->update(options);
        opts = options;
        
        if (not current_view)
            popAllViews();
            
        updateStates();
    }
    catch(const Conspire::Exception &e)
    {
        caughtException(e);
    }
}

void OptionsWidget::add(QString key)
{
    conspireDebug() << "OptionsWidget::add(" << key << ")";

    Options old_opts = opts;

    try
    {
        undo_stack->push( new OptionsAddCommand(this, key) );

        if (not current_view)
            popAllViews();
    }
    catch(const Conspire::Exception &e)
    {
        caughtException(e);
        setOptions(old_opts);
    }
}

void OptionsWidget::remove(QString key)
{
    Options old_opts = opts;

    try
    {
        undo_stack->push( new OptionsRemoveCommand(this, key) );

        if (not current_view)
            popAllViews();
    }
    catch(const Conspire::Exception &e)
    {
        caughtException(e);
        setOptions(old_opts);
    }
}

void OptionsWidget::update(QString key, Obj object)
{
    Options old_opts = opts;

    try
    {
        undo_stack->push( new OptionsUpdateCommand(this, key, object.asA<Value>()) );

        if (not current_view)
            popAllViews();
    }
    catch(const Conspire::Exception &e)
    {
        caughtException(e);
        setOptions(old_opts);
    }
}

void OptionsWidget::resizeEvent(QResizeEvent *event)
{
    if (event)
    {
        this->viewport()->resize(event->size());

        if (mainbar)
        {
            QPointF top_left = mapToScene(QRect(QPoint(0,0), QSize(2, 2)))
                                    .boundingRect().center();
                                    
            mainbar->setPos(top_left);
        }
    
        if (current_view)
        {
            current_view->resize(this->viewport()->size());

            setSceneRect( QRectF(0, 0, current_view->geometry().width(), 
                                       current_view->geometry().height()) );
        }
    }
}

void OptionsWidget::keyPressEvent(QKeyEvent *event)
{
    if (event)
    {
        if (event->key() == ::Qt::Key_Escape)
        {
            this->popView();
        }
        else
        {
            QGraphicsView::keyPressEvent(event);
        }
    }
}

/** Internal function used to update the states of the add, back, forward etc. */
void OptionsWidget::updateStates()
{
    OptionsEditView *e = dynamic_cast<OptionsEditView*>(current_view.data());
    
    if (e)
        emit( canAddChanged( e->options().canAddValues() ) );
    else
        emit( canAddChanged(false) );

    emit( canBackChanged( not view_history.isEmpty() ) );
    emit( canForwardChanged( not view_future.isEmpty() ) );
}

/** This function is called whenever any of the child views are deleted */
void OptionsWidget::viewDeleted(QObject *obj)
{
    conspireDebug() << "checkViews()";
    
    if (not current_view)
        this->popView();
}

/** Switch to a new view, pushing the old view onto the stack */
void OptionsWidget::pushView(QGraphicsWidget *v, bool clear_future)
{
    if (clear_future)
        view_future.clear();

    connect(v, SIGNAL(destroyed(QObject*)), this, SLOT(viewDeleted(QObject*)));

    if (v)
    {
        QParallelAnimationGroup *g = new QParallelAnimationGroup();
    
        if (current_view)
        {
            view_history.push(current_view);

            QPropertyAnimation *anim = new QPropertyAnimation(current_view,
                                                   "opacity");
                                                   
            anim->setDuration(1000);
            anim->setStartValue(1.0);
            anim->setEndValue(0.0);
            
            g->addAnimation(anim);
        }

        //make sure that the view has the right size
        v->setOpacity(0);
        v->setGeometry(0, 0, viewport()->width(), viewport()->height());

        //set up animating turning on the view
        {
            QPropertyAnimation *anim = new QPropertyAnimation(v, "opacity");
            
            anim->setDuration(1000);
            anim->setStartValue(0.0);
            anim->setEndValue(1.0);
            
            g->addAnimation(anim);
        }
        
        //now set up animating sliding to the new view
        {
            QPropertyAnimation *anim = new QPropertyAnimation(v, "x");
            
            anim->setDuration(1000);
            anim->setStartValue(2 * this->width());
            anim->setEndValue(0.0);
            anim->setEasingCurve(QEasingCurve::InOutBack);
            
            g->addAnimation(anim);
        }
        
        current_view = v;
        
        g->start(QAbstractAnimation::DeleteWhenStopped);

        updateStates();
    }
}

/** Switch back to the old view */
QGraphicsWidget* OptionsWidget::popView()
{
    if (view_history.isEmpty())
        return 0;

    QGraphicsWidget *old_view = 0;

    QParallelAnimationGroup *g = new QParallelAnimationGroup();

    if (current_view)
    {
        old_view = current_view;
        
        QPropertyAnimation *anim = new QPropertyAnimation(old_view, "opacity");
        anim->setDuration(1000);
        anim->setStartValue(1.0);
        anim->setEndValue(0.0);
        g->addAnimation(anim);
        
        {
            QPropertyAnimation *anim = new QPropertyAnimation(old_view, "x");
            anim->setDuration(1000);
            anim->setStartValue(0.0);
            anim->setEndValue(2 * this->width());
            anim->setEasingCurve(QEasingCurve::InOutBack);
            g->addAnimation(anim);
        }
    }
        
    current_view = view_history.pop();
    
    if (not current_view)
        current_view = top_view;
 
    current_view->setOpacity(0);
    current_view->setGeometry(0, 0, viewport()->width(), viewport()->height());
                      
    {
        QPropertyAnimation *anim = new QPropertyAnimation(current_view, "opacity");
        anim->setDuration(1000);
        anim->setStartValue(0.0);
        anim->setEndValue(1.0);
        g->addAnimation(anim);
    }

    g->start(QAbstractAnimation::DeleteWhenStopped);

    view_future.push(old_view);

    updateStates();

    return old_view;
}

/** Pop all views and return to the top-level view */
void OptionsWidget::popAllViews()
{
    if (current_view == top_view)
        return;

    QGraphicsWidget *old_view = 0;

    QParallelAnimationGroup *g = new QParallelAnimationGroup();

    if (current_view)
    {
        old_view = current_view;
        
        QPropertyAnimation *anim = new QPropertyAnimation(old_view, "opacity");
        anim->setDuration(1000);
        anim->setStartValue(1.0);
        anim->setEndValue(0.0);
        g->addAnimation(anim);
        
        {
            QPropertyAnimation *anim = new QPropertyAnimation(old_view, "x");
            anim->setDuration(1000);
            anim->setStartValue(0.0);
            anim->setEndValue(2 * this->width());
            anim->setEasingCurve(QEasingCurve::InOutBack);
            g->addAnimation(anim);
        }
    }
        
    current_view = top_view;
    
#warning NEED TO MOVE VIEW HISTORY INTO VIEW FUTURE
    view_history.clear();

    current_view->setOpacity(0);
    current_view->setGeometry(0, 0, viewport()->width(), viewport()->height());

    {
        QPropertyAnimation *anim = new QPropertyAnimation(current_view, "opacity");
        anim->setDuration(1000);
        anim->setStartValue(0.0);
        anim->setEndValue(1.0);
        g->addAnimation(anim);
    }

    updateStates();

    g->start(QAbstractAnimation::DeleteWhenStopped);
}

/** Internal function used to construct the widget */
void OptionsWidget::build()
{
    this->setMinimumSize(600,400);

    this->setScene( new QGraphicsScene(this) );
    
    top_view = new NewOptionsWidget();
    top_view->setGeometry(0, 0, viewport()->width(), viewport()->height());
    top_view->setZValue(0);
    current_view = top_view;
    
    setSceneRect( QRectF(0, 0, current_view->geometry().width(), 
                               current_view->geometry().height()) );

    this->scene()->addItem(top_view);

    undo_stack = new QUndoStack(this);

    mainbar = new MainBar();
    mainbar->setPos(0, 0);
    mainbar->setZValue(1);

    this->scene()->addItem(mainbar);

    connect(mainbar, SIGNAL(add()), this, SLOT(add()));
    connect(mainbar, SIGNAL(undo()), this, SLOT(undo()));
    connect(mainbar, SIGNAL(redo()), this, SLOT(redo()));
    connect(mainbar, SIGNAL(back()), this, SLOT(back()));
    connect(mainbar, SIGNAL(forward()), this, SLOT(forward()));
    
    connect(this, SIGNAL(canBackChanged(bool)), mainbar, SLOT(canBackChanged(bool)));
    connect(this, SIGNAL(canForwardChanged(bool)),
            mainbar, SLOT(canForwardChanged(bool)));
    connect(this, SIGNAL(canAddChanged(bool)), mainbar, SLOT(canAddChanged(bool)));
    
    connect(undo_stack, SIGNAL(canUndoChanged(bool)),
            mainbar, SLOT(canUndoChanged(bool)));
    connect(undo_stack, SIGNAL(canRedoChanged(bool)),
            mainbar, SLOT(canRedoChanged(bool)));
    
    connect(undo_stack, SIGNAL(undoTextChanged(QString)),
            mainbar, SLOT(undoTextChanged(QString)));
    connect(undo_stack, SIGNAL(redoTextChanged(QString)),
            mainbar, SLOT(redoTextChanged(QString)));
    
    connect(top_view, SIGNAL(add(QString)), this, SLOT(add(QString)));
    connect(top_view, SIGNAL(remove(QString)), this, SLOT(remove(QString)));
    connect(top_view, SIGNAL(update(QString,Obj)), this, SLOT(update(QString,Obj)));
    connect(top_view, SIGNAL(help(Option)), this, SLOT(help(Option)));
}

/** Called to undo the last action */
void OptionsWidget::undo()
{
    if (undo_stack->canUndo())
        undo_stack->undo();
}

/** Called to redo the last undone action */
void OptionsWidget::redo()
{
    if (undo_stack->canRedo())
        undo_stack->redo();
}

/** Called to return back to the previous screen */
void OptionsWidget::back()
{
    this->popView();
}

/** Called to move forward to the next screen */
void OptionsWidget::forward()
{
    if (view_future.isEmpty())
        return;
        
    QPointer<QGraphicsWidget> w = view_future.pop();
    
    if (not w)
    {
        //this widget has died - try to get another one
        this->forward();
        return;
    }
    
    this->pushView(w, false);
}

/** Function called when an exception is caught */
void OptionsWidget::caughtException(const Exception &e)
{
    QLabel *l = new QLabel();
    l->setWordWrap(true);
    l->setGeometry(0, 0, viewport()->width(), viewport()->height());
    l->setText(e.toString());
    
    QGraphicsProxyWidget *l_proxy = new QGraphicsProxyWidget();
    l_proxy->setWidget(l);
    l_proxy->setOpacity(0.0);
    
    scene()->addItem(l_proxy);
    
    // THIS LEAKS MEMORY !!!
    
    this->pushView(l_proxy);
}

/** Called to add a new option value to the current view */
void OptionsWidget::add()
{
    conspireDebug() << "ADD";
    AddWidget *a = new AddWidget(opts);
    
    connect(a, SIGNAL(add(QString)), this, SLOT(add(QString)));
    
    scene()->addItem(a);
    
    // The AddWidget will automatically delete itself
    // when it is no longer needed
    
    this->pushView(a);
}

/** Called to view the help for the passed option */
void OptionsWidget::help(Option option)
{
    conspireDebug() << "HELP" << option.description();
}
