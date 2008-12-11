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

#include "molpro.h"
#include "qmpotential.h"
#include "latticecharges.h"

#include "SireMol/element.h"

#include "SireBase/tempdir.h"
#include "SireBase/findexe.h"

#include "SireUnits/units.h"

#include "tostring.h"

#include "SireError/errors.h"

#include "SireStream/datastream.h"
#include "SireStream/shareddatastream.h"

#include <QDebug>

using namespace Squire;
using namespace SireFF;
using namespace SireMol;
using namespace SireVol;
using namespace SireBase;
using namespace SireUnits;
using namespace SireStream;
using namespace SireStream;

static const RegisterMetaType<Molpro> r_molpro;

/** Serialise to a binary datastream */
QDataStream SQUIRE_EXPORT &operator<<(QDataStream &ds, const Molpro &molpro)
{
    writeHeader(ds, r_molpro, 1);
    
    SharedDataStream sds(ds);
    
    sds << molpro.env_variables
        << molpro.molpro_exe << molpro.basis_set << molpro.qm_method
        << molpro.energy_template << molpro.force_template
        << molpro.total_charge << molpro.memory_requirement;
        
    return ds;
}

/** Extract from a binary datastream */
QDataStream SQUIRE_EXPORT &operator>>(QDataStream &ds, Molpro &molpro)
{
    VersionID v = readHeader(ds, r_molpro);
    
    if (v == 1)
    {
        SharedDataStream sds(ds);
        
        sds >> molpro.env_variables
            >> molpro.molpro_exe >> molpro.basis_set >> molpro.qm_method
            >> molpro.energy_template >> molpro.force_template
            >> molpro.total_charge >> molpro.memory_requirement;
    }
    else
        throw version_error(v, "1", r_molpro, CODELOC);
        
    return ds;
}

static const QString default_energy_template =
       "memory,@MEMORY_REQUIREMENT@\n"
       "geomtyp=xyz\n"
       "geometry={ NOSYM, NOORIENT,\n"
       "@NUM_QM_ATOMS@  ! number of atoms\n"
       "This is an auto-generated molpro command file generated using Sire\n"
       "@QM_COORDS@\n"
       "}\n\n"
       "set,CHARGE=@QM_CHARGE@\n"
       "lattice, NUCONLY\n"
       "BEGIN_DATA\n"
       "@LATTICE_POINTS@\n"
       "END\n\n"
       "basis=@BASIS_SET@\n\n"
       "@QM_METHOD@;SIRE_FINAL_ENERGY=energy\n";

static const QString default_force_template = "! NEEDS TO BE WRITTEN";

/** Constructor */
Molpro::Molpro() 
       : ConcreteProperty<Molpro,QMProgram>(),
         basis_set("vdz"), qm_method("HF"),
         energy_template(default_energy_template),
         force_template(default_force_template),
         total_charge(0), memory_requirement( 8 * 1024 * 1024 * 4 )
{}

/** Construct, passing in the location of the Molpro executable */
Molpro::Molpro(const QString &molpro)
       : ConcreteProperty<Molpro,QMProgram>(),
         basis_set("vdz"), qm_method("HF"),
         energy_template(default_energy_template),
         force_template(default_force_template),
         total_charge(0), memory_requirement( 8 * 1024 * 1024 * 4 )
{
    this->setExecutable(molpro);
}

/** Copy constructor */
Molpro::Molpro(const Molpro &other)
       : ConcreteProperty<Molpro,QMProgram>(other),
         env_variables(other.env_variables), molpro_exe(other.molpro_exe),
         basis_set(other.basis_set), qm_method(other.qm_method),
         energy_template(other.energy_template),
         force_template(other.force_template),
         total_charge(other.total_charge),
         memory_requirement(other.memory_requirement)
{}

/** Destructor */
Molpro::~Molpro()
{}

/** Copy assignment operator */
Molpro& Molpro::operator=(const Molpro &other)
{
    if (this != &other)
    {
        env_variables = other.env_variables;
        molpro_exe = other.molpro_exe;
        basis_set = other.basis_set;
        qm_method = other.qm_method;
        energy_template = other.energy_template;
        force_template = other.force_template;
        total_charge = other.total_charge;
        memory_requirement = other.memory_requirement;
    }
    
    return *this;
}

/** Comparison operator */
bool Molpro::operator==(const Molpro &other) const
{
    return this == &other or 
           (env_variables == other.env_variables and
            molpro_exe == other.molpro_exe and
            basis_set == other.basis_set and
            qm_method == other.qm_method and
            energy_template == other.energy_template and
            force_template == other.force_template and
            total_charge == other.total_charge and
            memory_requirement == other.memory_requirement);
}

/** Comparison operator */
bool Molpro::operator!=(const Molpro &other) const
{
    return not this->operator==(other);
}

/** Set the Molpro executable (full path and also arguments) to be used */
void Molpro::setExecutable(const QString &molpro_executable)
{
    molpro_exe = molpro_executable;
}

/** Return the executable (full path and also arguments) to be used. This
    is null if the executable is searched for in the path */
QString Molpro::executable() const
{
    return molpro_exe;
}

/** Set the environmental variable 'variable' to have the value 'value'
    when the Molpro executable is run. This replaces any existing
    value of this environmental variable */
void Molpro::setEnvironment(const QString &variable, const QString &value)
{
    env_variables[variable] = value;
}

/** Return all of the environmental variables that are to be set explicitly
    when Molpro is run. This does not include any environmental variables
    that have not been explicitly set, but do have values */
const QHash<QString,QString>& Molpro::environment() const
{
    return env_variables;
}

/** Return the value of the explicitly set environmental variable 'variable'.
    A null string is returned if this variable has not been set 
    explicitly (this does not mean the variable doesn't exist - merely
    that a specific value has not been set) */
QString Molpro::environment(const QString &variable) const
{
    return env_variables.value(variable);
}

/** Set the string representation of the basis set to be used
    during the calculation. This will substitute for @BASIS_SET@ in 
    the energy and force command file templates. */
void Molpro::setBasisSet(const QString &basis)
{
    basis_set = basis;
}

/** Return the basis set to be used during the calculation. This
    substitutes for @BASIS_SET@ in the energy and force command
    file templates */
const QString& Molpro::basisSet() const
{
    return basis_set;
}

/** Set the QM method to be used to calculate the energy or
    force (e.g. HF, df-ks,b3lyp). This will substitute for 
    @QM_METHOD@ in the energy and force command file templates */
void Molpro::setMethod(const QString &method)
{
    qm_method = method;
}

/** Return the QM method to be used to calculate the energy or
    force (e.g. HF, df-ks,b3lyp). This will substitute for 
    @QM_METHOD@ in the energy and force command file templates */
const QString& Molpro::method() const
{
    return qm_method;
}

/** Set the total charge of the system (in unit charges) */
void Molpro::setTotalCharge(int charge)
{
    total_charge = charge;
}

/** Return the total charge of the system */
int Molpro::totalCharge() const
{
    return total_charge;
}

/** Set the memory requirement (in bytes) that will be reserved for use
    by molpro. You will need to set this if the default amount
    (8000000 floating point words ~ 32MB) is not enough. Be careful
    not to allocate too much though or you will be swapping all of the time */
void Molpro::setMemoryRequirement(int nbytes)
{
    if (nbytes > 4)
        memory_requirement = nbytes;
}

/** Return the memory requirement that has been set for this job (in bytes) */
int Molpro::memoryRequirement() const
{
    return memory_requirement;
}

/** Set the template for the command file to be used to get
    Molpro to calculate an energy. The following tags will
    be substituted in the template;
    
    @BASIS_SET@          - the desired basis set
    @QM_METHOD@          - the desired QM method (e.g. HF)
    @QM_COORDS@          - the list of elements and coordinates of QM atoms
    @QM_CHARGE@          - the total charge of the system
    @NUM_QM_ATOMS@       - the number of QM atoms
    @LATTICE_POINTS@     - the coordinates and charges of the lattice points (MM atoms)
    @NUM_LATTICE_POINTS@ - the number of lattice points (MM atoms)
*/
void Molpro::setEnergyTemplate(const QString &cmd_template)
{
    energy_template = cmd_template;
}

/** Return the template for the command file to be used to get Molpro
    to calculate the energy. */
const QString& Molpro::energyTemplate() const
{
    return energy_template;
}

/** Set the template for the command file to be used to get
    Molpro to calculate the forces. The following tags will
    be substituted in the template;
    
    @BASIS_SET@          - the desired basis set
    @QM_METHOD@          - the desired QM method (e.g. HF)
    @QM_COORDS@          - the list of elements and coordinates of QM atoms
    @QM_CHARGE@          - the total charge of the system
    @NUM_QM_ATOMS@       - the number of QM atoms
    @LATTICE_POINTS@     - the coordinates and charges of the lattice points (MM atoms)
    @NUM_LATTICE_POINTS@ - the number of lattice points (MM atoms)
*/
void Molpro::setForceTemplate(const QString &cmd_template)
{
    force_template = cmd_template;
}

/** Return the template for the command file to be used to get Molpro
    to calculate the forces. */
const QString& Molpro::forceTemplate() const
{
    return force_template;
}

/** Function used to substitute in the atom and lattice coordinates
    into the provided molpro command template */
QString Molpro::createCommandFile(QString cmd_template,
                                  const QMPotential::Molecules &molecules,
                                  const LatticeCharges &lattice_charges) const
{
    //replace the easy things...
    cmd_template.replace( QLatin1String("@BASIS_SET@"), 
                          basis_set, Qt::CaseInsensitive );
    
    int required_memory_in_words = memory_requirement / 4;
    
    if (required_memory_in_words > 1000000)
    {
        cmd_template.replace( QLatin1String("@MEMORY_REQUIREMENT@"),
                              QString("%1,M").arg(required_memory_in_words / 1000000) );
    }
    else if (required_memory_in_words > 1000)
    {
        cmd_template.replace( QLatin1String("@MEMORY_REQUIREMENT@"),
                              QString("%1,K").arg(required_memory_in_words / 1000 ) );
    }
    else
    {
        cmd_template.replace( QLatin1String("@MEMORY_REQUIREMENT@"),
                              QString("%1").arg(required_memory_in_words) );
    }
    
    cmd_template.replace( QLatin1String("@QM_METHOD@"),
                          qm_method, Qt::CaseInsensitive );

    cmd_template.replace( QLatin1String("@QM_CHARGE@"),
                          QString::number(total_charge), Qt::CaseInsensitive );

    if (lattice_charges.isEmpty())
    {
        //there are no lattice charges
        cmd_template.replace( QLatin1String("@NUM_LATTICE_POINTS@"),
                              QString::number(0), Qt::CaseInsensitive );
                          
        cmd_template.replace( QLatin1String("@LATTICE_POINTS@"),
                              QLatin1String(" "), Qt::CaseInsensitive );
    }
                  
    //now build the list of all of the atoms
    QStringList atom_coords;
    
    int nmols = molecules.count();
    const QMPotential::Molecule *molecules_array = molecules.moleculesByIndex()
                                                            .constData();
                                                            
    for (int i=0; i<nmols; ++i)
    {
        const QMPotential::Molecule &molecule = molecules_array[i];
        
        //loop through the atoms...
        const CoordGroupArray &coords = molecule.coordinates();
        const PackedArray2D<Element> &elements = molecule.parameters().atomicParameters();
        
        int natoms = coords.nCoords();
        
        BOOST_ASSERT( natoms == elements.nValues() );
        
        const Vector *coords_array = coords.constCoordsData();
        const Element *elements_array = elements.constValueData();
        
        for (int j=0; j<natoms; ++j)
        {
            const Element &element = elements_array[j];
            
            if (element.nProtons() > 0)
            {
                //this is not a dummy atom!
                const Vector &c = coords_array[j];
    
                atom_coords.append( QString("%1,%2,%3,%4")
                                        .arg(element.symbol(),
                                             QString::number(c.x(), 'f', 8),
                                             QString::number(c.y(), 'f', 8),
                                             QString::number(c.z(), 'f', 8) ) );
            }
        }
    }
    
    cmd_template.replace( QLatin1String("@NUM_QM_ATOMS@"),
                          QString::number(atom_coords.count()), Qt::CaseInsensitive );
                          
    cmd_template.replace( QLatin1String("@QM_COORDS@"),
                          atom_coords.join("\n"), Qt::CaseInsensitive );
    
    //put the lattice charges in now (as they can make the command
    //file *very* large)
    if (not lattice_charges.isEmpty())
    {
        int ncharges = lattice_charges.nCharges();
        const LatticeCharge *charges_array = lattice_charges.constData();
        
        QStringList charges;
        
        for (int i=0; i<ncharges; ++i)
        {
            const LatticeCharge &charge = charges_array[i];
        
            if (charge.charge() != 0)
                charges.append( QString("%1,%2,%3,%4")
                                    .arg(QString::number(charge.x(), 'f', 8),
                                         QString::number(charge.y(), 'f', 8),
                                         QString::number(charge.z(), 'f', 8),
                                         QString::number(charge.charge(), 'f', 8) ) );
        }
        
        cmd_template.replace( QLatin1String("@NUM_LATTICE_POINTS@"),
                              QString::number(charges.count()), Qt::CaseInsensitive );
                          
        cmd_template.replace( QLatin1String("@LATTICE_POINTS@"),
                              charges.join("\n"), Qt::CaseInsensitive );
    }
                                       
    return cmd_template;
}

/** Return the command file that will be used to calculate the energy of the 
    molecules in 'molecules' */
QString Molpro::energyCommandFile(const QMPotential::Molecules &molecules) const
{
    return this->createCommandFile(energy_template, molecules);
}

/** Return the command file that will be used to calculate the energy of the 
    molecules in 'molecules' in the field of point charges in 'lattice_charges' */
QString Molpro::energyCommandFile(const QMPotential::Molecules &molecules,
                                  const LatticeCharges &lattice_charges) const
{
    return this->createCommandFile(energy_template, molecules, lattice_charges);
}

/** Return the command files that will be used to calculate the forces on the  
    atoms of the molecules in 'molecules' */
QString Molpro::forceCommandFile(const QMPotential::Molecules &molecules) const
{
    return this->createCommandFile(force_template, molecules);
}

/** Return the command files that will be used to calculate the forces on the  
    atoms of the molecules in 'molecules' in the field of point charges
    in 'lattice_charges' - this also calculates the forces on those
    point charges */
QString Molpro::forceCommandFile(const QMPotential::Molecules &molecules,
                                 const LatticeCharges &lattice_charges) const
{
    return this->createCommandFile(force_template, molecules, lattice_charges);
}

/** Extract the energy from the molpro output in 'molpro_output' */
double Molpro::extractEnergy(QFile &molpro_output) const
{
    QTextStream ts(&molpro_output);

    QRegExp regexp("SIRE_FINAL_ENERGY\\s*=\\s*([-\\d\\.]+)");

    QStringList lines;
    
    while (not ts.atEnd())
    {
        QString line = ts.readLine();
        lines.append(line);

        if (regexp.indexIn(line) != -1)
        {
            //we've found the SIRE_FINAL_ENERGY line
            QString num = regexp.cap(1);
        
            bool ok;
        
            double nrg = num.toDouble(&ok);
        
            if (not ok)
                throw SireError::process_error( QObject::tr(
                    "The energy obtained from Molpro is garbled (%1) - %2.")
                        .arg(regexp.cap(1), regexp.cap(0)), CODELOC );
        
            //the energy is in hartrees - convert it to kcal per mol
            return nrg * hartree;
        }
    }

    //ok - something went wrong - we couldn't find the output
    // Lets try reading the file again (in case it wasn't fully written
    // to the disc)
    QTextStream ts2(&molpro_output);

    lines = QStringList();
    
    while (not ts2.atEnd())
    {
        QString line = ts2.readLine();
        lines.append(line);

        if (regexp.indexIn(line) != -1)
        {
            //we've found the SIRE_FINAL_ENERGY line
            QString num = regexp.cap(1);
        
            bool ok;
        
            double nrg = num.toDouble(&ok);
        
            if (not ok)
                throw SireError::process_error( QObject::tr(
                    "The energy obtained from Molpro is garbled (%1) - %2.")
                        .arg(regexp.cap(1), regexp.cap(0)), CODELOC );
        
            //the energy is in hartrees - convert it to kcal per mol
            return nrg * hartree;
        }
    }

    //the SIRE_FINAL_ENERGY line was not found!
    throw SireError::process_error( QObject::tr(
            "Could not find the total energy in the molpro output!\n"
                "%1").arg(lines.join("\n")), CODELOC );
}

/** Internal function used to write the shell script that is used to
    run the molpro job and collect the output
*/
QString Molpro::writeShellFile(const TempDir &tempdir) const
{
    QString cmdfile = QString("%1/run_molpro.cmd").arg(tempdir.path());
    
    QFile f(cmdfile);
    f.open(QIODevice::WriteOnly);
    
    QTextStream ts(&f);

    //set the environmental variables of the job
    for (QHash<QString,QString>::const_iterator it = env_variables.constBegin();
         it != env_variables.constEnd();
         ++it)
    {
        ts << it.key() << "=\"" << it.value() << "\"\n";
    }

    //set the script to change into the run directory of the job
    ts << QString("\ncd %1").arg(tempdir.path()) << "\n\n";

    //write the line used to run molpro
    if (molpro_exe.isEmpty())
    {
        //the user hasn't specified a molpro executable - try to find one
        QString found_molpro = SireBase::findExe("molpro").absoluteFilePath();
        ts << QString("%1 -d %2 < molpro_input > molpro_output\n")
                    .arg(found_molpro, tempdir.path());
    }
    else
        ts << QString("%1 -d %2 < molpro_input > molpro_output\n")
                        .arg(molpro_exe, tempdir.path());
    
    f.close();
    
    return cmdfile;
}

static QByteArray readAll(const QString &file)
{
    QFile f(file);
    
    if (f.open(QIODevice::ReadOnly))
    {
        return f.readAll();
    }
    else
        return QByteArray();
}

/** Return the energy calculate according to the Molpro command
    file 'cmd_file' (this is the contents of the file, not
    the path to the file) */
double Molpro::calculateEnergy(const QString &cmdfile, int ntries) const
{
    //create a temporary directory in which to run Molpro
    QString tmppath = env_variables.value("TMPDIR");
    
    if (tmppath.isEmpty())
        tmppath = QDir::temp().absolutePath();

    TempDir tmpdir(tmppath);

    //write the file processed by the shell used to run the job
    QString shellfile = this->writeShellFile(tmpdir);

    {
        QFile f( QString("%1/molpro_input").arg(tmpdir.path()) );
        f.open( QIODevice::WriteOnly );
   
        //write the command file
        f.write( cmdfile.toAscii() );
        f.close();
    }

    //now run the command
    if (std::system(0))
    {
        std::system( QString("source %1").arg(shellfile).toAscii().constData() );
    }
    else
        throw SireError::unsupported( QObject::tr(
            "The operating system does not support the use of std::system(). "
            "It is not possible to use Molpro to calculate QM energies."),
                CODELOC );

    //read all of the output
    QFile f( QString("%1/molpro_output").arg(tmpdir.path()) );
    
    if ( not (f.exists() and f.open(QIODevice::ReadOnly)) )
    {
        QByteArray shellcontents = ::readAll(shellfile);
        QByteArray cmdcontents = ::readAll(QString("%1/molpro_input").arg(tmpdir.path()));
    
        throw SireError::process_error( QObject::tr(
            "There was an error running the Molpro job - no output was created.\n"
            "The shell script used to run the job was;\n"
            "*****************************************\n"
            "%1\n"
            "*****************************************\n"
            "The molpro input used to run the job was;\n"
            "*****************************************\n"
            "%2\n"
            "*****************************************\n")
                .arg( QLatin1String(shellcontents),
                      QLatin1String(cmdcontents) ), CODELOC );
    }


    try
    {
        //parse the output to get the energy
        return this->extractEnergy(f);
    }
    catch(...)
    {
        qDebug() << "Molpro process error. Number of remaining attempts = " << ntries;
    
        if (ntries <= 0)
            //don't bother trying again - it's not going to work!
            throw;
            
        //give it one more go - you never know, it may work
        return this->calculateEnergy(cmdfile, ntries-1);
    }
}

/** Run Molpro and use it to calculate the energy of the molecules in 
    'molecules'. This blocks until Molpro has completed */
double Molpro::calculateEnergy(const QMPotential::Molecules &molecules,
                               int ntries) const
{
    //create the command file to be used by Molpro
    QString cmdfile = this->energyCommandFile(molecules);
    
    return this->calculateEnergy(cmdfile, ntries);
}

/** Calculate the Molpro QM energy of the molecules in 'molecules'
    in the field of point charges in 'lattice_charges' */
double Molpro::calculateEnergy(const QMPotential::Molecules &molecules,
                               const LatticeCharges &lattice_charges,
                               int ntries) const
{
    //create the command file to be used by Molpro
    QString cmdfile = this->energyCommandFile(molecules, lattice_charges);
    
    return this->calculateEnergy(cmdfile, ntries);
}
