/********************************************\
  *
  *  Sire - Molecular Simulation Framework
  *
  *  Copyright (C) 2008  Christopher Woods
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

#ifndef SQUIRE_QMMMELECEMBEDPOTENTIAL_H
#define SQUIRE_QMMMELECEMBEDPOTENTIAL_H

#include "qmmmpotential.h"

#include "qmpotential.h"

#include "SireMM/intercoulombff.h"

SIRE_BEGIN_HEADER

namespace Squire
{
class QMMMElecEmbedPotential;
}

QDataStream& operator<<(QDataStream&, const Squire::QMMMElecEmbedPotential&);
QDataStream& operator>>(QDataStream&, Squire::QMMMElecEmbedPotential&);

namespace Squire
{

/** This is a QM/MM potential that uses electrostatic embedding to 
    allow the MM point charges to polarise the QM wavefunction
    
    @author Christopher Woods
*/
class SQUIRE_EXPORT QMMMElecEmbedPotential
           : public QMMMPotential<QMPotential, SireMM::InterCoulombPotential>
{

friend QDataStream& ::operator<<(QDataStream&, const QMMMElecEmbedPotential&);
friend QDataStream& ::operator>>(QDataStream&, QMMMElecEmbedPotential&);

public:
    typedef QMEnergy Energy;
    typedef Energy::Components Components;

    QMMMElecEmbedPotential();
    QMMMElecEmbedPotential(const QMMMElecEmbedPotential &other);
    
    ~QMMMElecEmbedPotential();
    
    QMMMElecEmbedPotential& operator=(const QMMMElecEmbedPotential &other);
    
    static const char* typeName()
    {
        return "Squire::QMMMElecEmbedPotential";
    }
    
    const char* what() const
    {
        return QMMMElecEmbedPotential::typeName();
    }

    bool setProperty(const QString &name, const Property &property);
    const Property& property(const QString &name) const;
    bool containsProperty(const QString &name) const;
    const Properties& properties() const;

    QString energyCommandFile(const QMMolecules &qmmols,
                              const MMMolecules &mmmols) const;
                              
    QString forceCommandFile(const QMMolecules &qmmols,
                             const MMMolecules &mmmols) const;
};

}

SIRE_END_HEADER

#endif
