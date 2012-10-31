#ifndef CONSPIRE_WORKVIEW_H
#define CONSPIRE_WORKVIEW_H
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

#include "Conspire/GUI/page.h"

#include "Conspire/option.h"

CONSPIRE_BEGIN_HEADER

#include <QGraphicsView>

//class QGraphicsView;

namespace Conspire
{
    /** This is a subclassing of QGraphicsView to give a view of the work
      * whose events can be grabbed and put to work using context sensitive
      * modes. **/
    class CONSPIRE_EXPORT WorkView : public QGraphicsView
    {
        Q_OBJECT
        
    public:
        WorkView();
        ~WorkView();

    protected slots:
        
    protected:
        void drawBackground(QPainter *painter, const QRectF &rect);
        void drawForeground(QPainter *painter, const QRectF &rect);
        void scaleView(qreal scaleFactor);
    };

}

CONSPIRE_END_HEADER

#endif
