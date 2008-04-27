/********************************************\
  *
  *  Sire - Molecular Simulation Framework
  *
  *  Copyright (C) 2007  Christopher Woods
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

#include "iobase.h"

#include "SireMol/mover.hpp"

#include "SireError/errors.h"

#include <QFile>

using namespace SireIO;
using namespace SireMol;
using namespace SireBase;

/** Constructor */
IOBase::IOBase() : PropertyBase()
{}

/** Copy constructor */
IOBase::IOBase(const IOBase &other) : PropertyBase(other)
{}

/** Destructor */
IOBase::~IOBase()
{}

/** Read all of the molecules contained in the file 'filename', using
    the (optional) passed properties in 'map', and returning a molecule
    group containing the molecules in the same order as they appear
    in the file, and with the molecule group name being 'filename' */
MoleculeGroup IOBase::read(const QString &filename, const PropertyMap &map) const
{
    //open the file for reading
    QFile fle(filename);
    if (not fle.open(QIODevice::ReadOnly|QIODevice::Text))
    {
        throw SireError::file_error(fle, CODELOC);
    }

    //get all of the data out of the file
    QByteArray contents = fle.readAll();
    fle.close();

    if (contents.isEmpty())
        return MoleculeGroup();

    MoleculeGroup molecules = this->read(contents, map);
    
    molecules.edit().setName(filename);

    return molecules;
}

/** Read all of the molecules contained in the IO device 'dev', using
    the (optional) passed properties in 'map', and returning a molecule
    group containing the molecules in the same order as they appear
    on the device, and with the molecule group name being taken
    from the device */
MoleculeGroup IOBase::read(QIODevice &dev, const PropertyMap &map) const
{
    if (not dev.isReadable())
    {
        throw SireError::io_error(QObject::tr("Cannot read molecules from an unreadble "
                                              "device!"), CODELOC );
    }

    //get all of the data from the device
    QByteArray contents = dev.readAll();

    if (contents.isEmpty())
        return MoleculeGroup();

    MoleculeGroup molecules = this->readMols(contents, map);

    //maybe one day try to name the group...

    return molecules;
}

/** Write the molecules in the passed group to the file called 'filename'.
    This writes the molecules in the same order as they appear in the
    passed group. */
void IOBase::write(const MolGroup &molgroup, const QString &filename,
                   const PropertyMap &map) const
{
    //write the group to a binary blob
    QByteArray data = this->writeMols(molgroup, map);

    //open a file into which to write the molecules
    QFile fle(filename);
    
    if (not fle.open(QIODevice::WriteOnly|QIODevice::Text))
    {
        throw SireError::file_error(fle);
    }
    
    //dump the blob into the file
    if (fle.write(data) == -1)
    {
        throw SireError::file_error(fle);
    }
}

/** Write the molecules in the passed group to the file called 'filename'. */
void IOBase::write(const Molecules &molecules, const QString &filename,
                   const PropertyMap &map) const
{
    //write the group to a binary blob
    QByteArray data = this->writeMols(molecules, map);

    //open a file into which to write the molecules
    QFile fle(filename);
    
    if (not fle.open(QIODevice::WriteOnly|QIODevice::Text))
    {
        throw SireError::file_error(fle);
    }
    
    //dump the blob into the file
    if (fle.write(data) == -1)
    {
        throw SireError::file_error(fle);
    }
}

/** Write the molecules in the passed group to the IO device 'dev'.
    This writes the molecules in the same order as they appear in the
    passed group. */
void IOBase::write(const MolGroup &molgroup, QIODevice &dev,
                   const PropertyMap &map) const
{
    //write the group to a binary blob
    QByteArray data = this->writeMols(molgroup, map);

    //dump the blob to the device
    if (dev.write(data) == -1)
    {
        throw SireError::file_error(dev.errorString(), CODELOC);
    }
}

/** Write the molecules in the passed group to the IO device 'dev'. */
void IOBase::write(const Molecules &molecules, QIODevice &dev,
                   const PropertyMap &map) const
{
    //write the group to a binary blob
    QByteArray data = this->writeMols(molecules, map);

    //dump the blob into the file
    if (dev.write(data) == -1)
    {
        throw SireError::file_error(dev.errorString(), CODELOC);
    }
}
