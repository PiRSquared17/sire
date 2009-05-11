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

#include <QFile>
#include <QDir>
#include <QByteArray>
#include <QTemporaryFile>

#include "trajectorymonitor.h"
#include "pdb.h"

#include "SireSystem/system.h"

#include "SireMol/moleculegroup.h"

#include "SireStream/datastream.h"
#include "SireStream/shareddatastream.h"

#include "SireError/errors.h"

#include <QDebug>

using namespace SireIO;
using namespace SireSystem;
using namespace SireMol;
using namespace SireBase;
using namespace SireStream;

static const RegisterMetaType<TrajectoryMonitor> r_trajmon;

static QByteArray readFromDisk(const QTemporaryFile &tmpfile)
{
    //open the packed data file using a separate file handle
    QFile f( tmpfile.fileName() );
    
    if (not f.open(QIODevice::ReadOnly | QIODevice::Unbuffered))
    {
        throw SireError::file_error(f, CODELOC);
    }
    
    //now read and return all of the data
    return f.readAll();
}

static QString getUserName()
{
    #ifdef Q_OS_UNIX
        return std::getenv("USER");
    #else
        return "USER";
    #endif
}

static boost::shared_ptr<QTemporaryFile> writeToDisk(const QByteArray &data,
                                                     const QString &temp_dir)
{
    QString save_dir;

    if (temp_dir.isEmpty())
    {
        //save to QDir::temp()
        save_dir = QDir::tempPath();
    }
    else
    {
        QDir tempdir(temp_dir);
    
        if (tempdir.exists())
            save_dir = tempdir.absolutePath();
        else
            save_dir = QDir::tempPath();
    }
    
    boost::shared_ptr<QTemporaryFile> tmp;
    
    tmp.reset( new QTemporaryFile( QString("%1/%2_sire_trajmonitor_XXXXXX.data")
                            .arg(save_dir, getUserName()) ) );
                            
    if (not tmp->open())
    {
        throw SireError::file_error( QObject::tr(
            "Cannot save the trajectory to disk as the temporary directory "
            "(%1) is not writable. Check that you have permission and "
            "there is sufficient space.")
                .arg(save_dir), CODELOC );
    }

    qint64 nbytes = tmp->write(data);

    if (nbytes == -1)
    {
        throw SireError::file_error( QObject::tr(
            "There was an error writing the trajectory to disk. "
            "Maybe you don't have write permission or maybe there is not "
            "enough disk space (we require at least %1 MB in %2)")
                .arg( data.count() / (1024.0*1024.0) )
                .arg(save_dir), CODELOC );
    }
    else if (nbytes != data.count())
    {
        throw SireError::file_error( QObject::tr(
            "There was a problem writing the trajectory to disk. "
            "Only %1 bytes were written, out of the necessary %2 bytes. "
            "Maybe the disk became full while it was being written (to %3)?")
                .arg(nbytes).arg(data.count())
                .arg(save_dir), CODELOC );
    }
    
    return tmp;
}

/** Serialise to a binary datastream */
QDataStream SIREIO_EXPORT &operator<<(QDataStream &ds, 
                                      const TrajectoryMonitor &trajmon)
{
    writeHeader(ds, r_trajmon, 1);
    
    SharedDataStream sds(ds);
    
    sds << trajmon.io_writer << trajmon.molgroup
        << trajmon.mol_properties << trajmon.temp_dir;

    //now write all of the frames
    sds << quint32( trajmon.traj_frames.count() );
    
    for (QList< boost::shared_ptr<QTemporaryFile> >::const_iterator it 
                                                   = trajmon.traj_frames.constBegin();
         it != trajmon.traj_frames.constEnd();
         ++it)
    {
        QByteArray data = ::readFromDisk( *(*it) );

        //don't shared-stream as this data cannot be shared, and 
        //pointer re-use may lead to fake shared streaming
        ds << data;
    }
        
    sds << static_cast<const SystemMonitor&>(trajmon);
        
    return ds;
}

/** Extract from a binary datastream */
QDataStream SIREIO_EXPORT &operator>>(QDataStream &ds, TrajectoryMonitor &trajmon)
{
    VersionID v = readHeader(ds, r_trajmon);
    
    if (v == 1)
    {
        SharedDataStream sds(ds);
        
        TrajectoryMonitor new_monitor;
        
        sds >> new_monitor.io_writer >> new_monitor.molgroup
            >> new_monitor.mol_properties >> new_monitor.temp_dir;

        //how many frames need to be read?
        quint32 nframes;
        sds >> nframes;
        
        for (quint32 i=0; i<nframes; ++i)
        {
            QByteArray data;
            
            //no shared streaming here
            ds >> data;
            
            new_monitor.traj_frames.append( ::writeToDisk(data, new_monitor.temp_dir) );
        }
            
        sds >> static_cast<SystemMonitor&>(trajmon);
    }
    else
        throw version_error(v, "1", r_trajmon, CODELOC);
        
    return ds;
}

/** Null constructor */
TrajectoryMonitor::TrajectoryMonitor()
                  : ConcreteProperty<TrajectoryMonitor,SystemMonitor>()
{}

/** Construct a monitor that monitors the trajectory of the molecules  
    in the molecule group 'molgroup'. This writes the trajectory using
    the PDB writer, and uses the (optionally) supplied property map to
    control what is written */
TrajectoryMonitor::TrajectoryMonitor(const MoleculeGroup &mol_group,
                                     const PropertyMap &map)
                  : ConcreteProperty<TrajectoryMonitor,SystemMonitor>(),
                    io_writer( PDB() ), molgroup(mol_group), mol_properties(map)
{}

/** Construct a monitor that monitors the trajectory of the molecules in 
    the molecule group 'molgroup', writing the trajectory using the 
    molecule write 'writer', and using the (optionally) supplied property
    map to control what is written */
TrajectoryMonitor::TrajectoryMonitor(const MoleculeGroup &mol_group, 
                                     const IOBase &writer,
                                     const PropertyMap &map)
                  : ConcreteProperty<TrajectoryMonitor,SystemMonitor>(),
                    io_writer(writer), molgroup(mol_group), mol_properties(map)
{}

/** Copy constructor */
TrajectoryMonitor::TrajectoryMonitor(const TrajectoryMonitor &other)
                  : ConcreteProperty<TrajectoryMonitor,SystemMonitor>(other),
                    io_writer(other.io_writer), molgroup(other.molgroup),
                    traj_frames(other.traj_frames), mol_properties(other.mol_properties),
                    temp_dir(other.temp_dir)
{}

/** Destructor */
TrajectoryMonitor::~TrajectoryMonitor()
{}

/** Copy assignment operator */
TrajectoryMonitor& TrajectoryMonitor::operator=(const TrajectoryMonitor &other)
{
    SystemMonitor::operator=(other);
    
    io_writer = other.io_writer;
    molgroup = other.molgroup;
    traj_frames = other.traj_frames;
    mol_properties = other.mol_properties;
    temp_dir = other.temp_dir;
    
    return *this;
}

/** Comparison operator */
bool TrajectoryMonitor::operator==(const TrajectoryMonitor &other) const
{
    return (this == &other) or
           (io_writer == other.io_writer and
            molgroup == other.molgroup and
            mol_properties == other.mol_properties and
            temp_dir == other.temp_dir and
            traj_frames == other.traj_frames and
            SystemMonitor::operator==(other));
}

/** Comparison operator */
bool TrajectoryMonitor::operator!=(const TrajectoryMonitor &other) const
{
    return not this->operator==(other);
}
    
/** Set the temporary directory used to store the trajectory as it 
    is being monitored during the simulation */
void TrajectoryMonitor::setTempDir(const QString &tempdir)
{
    temp_dir = tempdir;
}

static QString getFrameNumber(int i, int n)
{
    if (n == 1)
        return "";
    
    else
    {
        return QString("%1").arg( i, n, int(10), QLatin1Char('0'));
    }
}

static int nColumns(int n)
{
    if (n <= 0)
        return 0;
    else
        return int( log10(n) ) + 1;
}

/** Write the trajectory to disk, using the file template 'file_template'. 
    This writes each frame to a separate file, numbering them sequentially
    from 0, replacing "XXXXXX" with the frame number */
void TrajectoryMonitor::writeToDisk(const QString &file_template) const
{
    int nframes = traj_frames.count();
    int i = 0;
    int n = nColumns(nframes);

    for (QList< boost::shared_ptr<QTemporaryFile> >::const_iterator
                            it = traj_frames.constBegin();
         it != traj_frames.constEnd();
         ++it)
    {
        //read the file
        QByteArray data = ::readFromDisk( *(*it) );
        data = qUncompress(data);
    
        QString filename = file_template;
        
        if (filename.contains("XXXXXX"))
            filename.replace("XXXXXX", getFrameNumber(i, n));
        else
            filename += (getFrameNumber(i,n) + ".pdb");
            
        ++i;
        
        QFile f(filename);
        
        if (not f.open( QIODevice::WriteOnly ))
            throw SireError::file_error(f);
            
        qint64 nbytes = f.write(data);

        if (nbytes == -1)
        {
            throw SireError::file_error( QObject::tr(
                "There was an error writing a trajectory frame to disk. "
                "Maybe you don't have write permission or maybe there is not "
                "enough disk space (we require at least %1 MB for %2)")
                .arg(data.count() / (1024.0*1024.0) )
                .arg(filename), CODELOC );
        }
        else if (nbytes != data.count())
        {
            throw SireError::file_error( QObject::tr(
                "There was a problem writing the trajectory frame to disk. "
                "Only %1 bytes were written, out of the necessary %2 bytes. "
                "Maybe the disk became full while it was being written (to %3)?")
                    .arg(nbytes).arg(data.count())
                    .arg(filename), CODELOC );
        }
    }
}

/** Clear all statistics - this will clear all frames of the trajectory */
void TrajectoryMonitor::clearStatistics()
{
    traj_frames.clear();
}

/** Monitor the system, writing an additional frame of the trajectory
    to this monitor */
void TrajectoryMonitor::monitor(System &system)
{
    if (io_writer.isNull())
        //there is nothing to write
        return;

    if (system.contains( molgroup.number() ))
    {
        const MoleculeGroup &new_group = system[molgroup.number()];
        
        if (traj_frames.isEmpty() or new_group.version() != molgroup.version())
        {
            //we need to write a new frame
            QByteArray frame_data = io_writer->write(new_group, mol_properties);
            
            //compress the data and save it to a temporary file
            frame_data = qCompress(frame_data);
            
            //now save this data to a temporary file
            traj_frames.append( ::writeToDisk(frame_data, temp_dir) );
            
            //save the new state of the molecule group
            molgroup = new_group;
        }
        else
        {
            //copy the last frame
            traj_frames.append( traj_frames.last() );
        }
    }
}