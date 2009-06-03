/********************************************\
  *
  *  Sire - Molecular Simulation Framework
  *
  *  Copyright (C) 2009  Christopher Woods
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

#include <QMutex>

#include "gldisplaylist.h"

#include <QDebug>

namespace Spier
{
namespace detail
{

class GLDisplayListData
{
public:
    GLDisplayListData(const QGLContext *context);
    ~GLDisplayListData();
    
    void aboutToDelete(const QGLContext *context);
    
    GLuint list_id;
    const QGLContext *render_context;
    
    QMutex list_mutex;
    
    bool is_recording;
};

} // end of namespace detail
} // end of namespace Spier


using namespace Spier;
using namespace Spier::detail;

////////
//////// Implementation of GLDisplayListData
////////

GLDisplayListData::GLDisplayListData(const QGLContext *context)
                  : list_id(0), 
                    render_context(context),
                    list_mutex(QMutex::Recursive),
                    is_recording(false)
{
    if (context)
    {
        list_id = glGenLists(1);
                
        if (list_id)
        {
            is_recording = true;
            glNewList(list_id, GL_COMPILE_AND_EXECUTE);
            list_mutex.lock();
        }
    }
}

GLDisplayListData::~GLDisplayListData()
{
    if (is_recording)
        qDebug() << "Trying to delete a GLDisplayListData object that "
                 << "is being recorded!!!";

    if (list_id)
        glDeleteLists(list_id, 1);
}

void GLDisplayListData::aboutToDelete(const QGLContext *context)
{
    if (not list_id)
        return;

    if (render_context != context)
        qDebug() << "Uh oh - we are going to delete a display list in the"
                    "wrong context...";
                    
    glDeleteLists(list_id, 1);
    list_id = 0;
}

////////
//////// Implementation of GLDisplayRecorder
////////

/** Constructor */
GLDisplayListRecorder::GLDisplayListRecorder() : boost::noncopyable()
{}

/** Construct for the passed GLDisplayList */
GLDisplayListRecorder::GLDisplayListRecorder(GLDisplayList &list, 
                                             const QGLContext *context)
                      : boost::noncopyable()
{
    list.startRecording(context);
    recording_list = list;
}

/** Destructor */
GLDisplayListRecorder::~GLDisplayListRecorder()
{
    recording_list.stopRecording();
}

////////
//////// Implementation of GLDisplayList
////////

/** Constructor */
GLDisplayList::GLDisplayList()
{}

/** Copy constructor */
GLDisplayList::GLDisplayList(const GLDisplayList &other) : d(other.d)
{}

/** Destructor */
GLDisplayList::~GLDisplayList()
{}

/** Copy assignment operator */
GLDisplayList& GLDisplayList::operator=(const GLDisplayList &other)
{
    d = other.d;
    return *this;
}

/** Comparison operator */
bool GLDisplayList::operator==(const GLDisplayList &other) const
{
    return d.get() == other.d.get();
}

/** Comparison operator */
bool GLDisplayList::operator!=(const GLDisplayList &other) const
{
    return d.get() != other.d.get();
}

/** Play this display list */
bool GLDisplayList::play(const QGLContext *playing_context)
{
    if (d.get() != 0)
    {
        QMutexLocker lkr( &(d->list_mutex) );

        if (d->render_context == playing_context)
        {
            glCallList( d->list_id );
            return true;
        }
    }

    return false;
}

/** Start recording a new display list */
void GLDisplayList::startRecording(const QGLContext *render_context)
{
    if (d.unique())
    {
        d->aboutToDelete(render_context);
    }
    
    d.reset( new GLDisplayListData(render_context) );
}

/** Stop recording a new display list */
void GLDisplayList::stopRecording()
{
    if (d.get() != 0)
    {
        QMutexLocker lkr( &(d->list_mutex) );
        
        if (d->is_recording)
        {
            glEndList();
            d->is_recording = false;
            d->list_mutex.unlock();
        }
    }
}

/** Is this an empty list? */
bool GLDisplayList::isEmpty() const
{
    return d.get() == 0;
}

/** Return whether or not this is a valid display list for 
    the passed QGLContext */
bool GLDisplayList::isValidFor(const QGLContext *render_context) const
{
    return (d->render_context == render_context);
}

/** Return the render context for this display list */
const QGLContext* GLDisplayList::renderContext() const
{
    return d->render_context;
}

////////
//////// Implementation of GLDisplayListRegistry
////////

GLDisplayListRegistry::GLDisplayListRegistry()
{}

GLDisplayListRegistry::~GLDisplayListRegistry()
{
    //delete each display list
    display_lists.clear();
    lists_to_delete.clear();
}

Q_GLOBAL_STATIC( QMutex, dlistRegistryMutex );

boost::shared_ptr< QHash<const QGLContext*,GLDisplayListRegistry> > 
        GLDisplayListRegistry::registries;

/** Return the display list for the object with key 'key' in the 
    render context 'render_context' */
GLDisplayList GLDisplayListRegistry::getDisplayList(const QGLContext *render_context,
                                                    const void *key)
{
    QMutexLocker lkr( dlistRegistryMutex() );
    
    if (registries.get() == 0)
        registries.reset( new QHash<const QGLContext*,GLDisplayListRegistry>() );


    GLDisplayListRegistry &registry = (*registries)[render_context];

    registry.lists_to_delete.clear();
    return registry.display_lists.value(key, GLDisplayList());
}

/** Save the display list 'display_list' for the QGLContext render context
    'render_context' for the object with key 'key' */
void GLDisplayListRegistry::saveDisplayList(const QGLContext *render_context,
                                            const void *key,
                                            const GLDisplayList &display_list)
{
    if (display_list.renderContext() != render_context)
        return;

    QMutexLocker lkr( dlistRegistryMutex() );
    
    if (registries.get() == 0)
        registries.reset( new QHash<const QGLContext*,GLDisplayListRegistry>() );

    GLDisplayListRegistry &registry = (*registries)[render_context];

    registry.display_lists.insert(key, display_list);
    registry.lists_to_delete.clear();
}

void GLDisplayListRegistry::clearLists(const void *key)
{
    if (display_lists.contains(key))
    {
        lists_to_delete.append( display_lists.take(key) );
    }
}

/** Clear all display lists associated with the key 'key' */
void GLDisplayListRegistry::clearDisplayLists(const void *key)
{
    QMutexLocker lkr( dlistRegistryMutex() );
    
    if (registries.get() == 0)
        registries.reset( new QHash<const QGLContext*,GLDisplayListRegistry>() );

    for (QHash<const QGLContext*,GLDisplayListRegistry>::iterator
                                                it = registries->begin();
         it != registries->end();
         ++it)
    {
        it->clearLists(key);
    }
}

/** Clear all display lists associated with the render context 'render_context' */
void GLDisplayListRegistry::clearDisplayLists(const QGLContext *render_context)
{
    QMutexLocker lkr( dlistRegistryMutex() );
    
    if (registries.get() == 0)
        registries.reset( new QHash<const QGLContext*,GLDisplayListRegistry>() );

    registries->remove(render_context);
}

void GLDisplayListRegistry::copyList(const void *old_key, const void *new_key)
{
    if (display_lists.contains(old_key))
    {
        if (display_lists.contains(new_key))
            lists_to_delete.append( display_lists.take(new_key) );
    
        display_lists.insert( new_key, display_lists.value(old_key) );
    }
}

/** Copy the display lists of key 'old_key' to the display list 'new_key' */
void GLDisplayListRegistry::copyDisplayList(const void *old_key, const void *new_key)
{
    if (old_key == new_key)
        return;

    QMutexLocker lkr( dlistRegistryMutex() );
    
    if (registries.get() == 0)
        registries.reset( new QHash<const QGLContext*,GLDisplayListRegistry>() );

    for (QHash<const QGLContext*,GLDisplayListRegistry>::iterator
                                                it = registries->begin();
         it != registries->end();
         ++it)
    {
        it->copyList(old_key, new_key);
    }
}