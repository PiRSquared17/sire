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

#include <cmath>

#include "montecarlo.h"

#include "SireFF/forcefields.h"

#include "SireUnits/units.h"
#include "SireUnits/temperature.h"

#include "SireStream/datastream.h"
#include "SireStream/shareddatastream.h"

#include <QDebug>

using namespace SireMove;
using namespace SireUnits;
using namespace SireSystem;
using namespace SireFF;
using namespace SireUnits::Dimension;
using namespace SireStream;

static const RegisterMetaType<MonteCarlo> r_mc(MAGIC_ONLY, "SireMove::MonteCarlo");

/** Serialise to a binary data stream */
QDataStream SIREMOVE_EXPORT &operator<<(QDataStream &ds, const MonteCarlo &mc)
{
    writeHeader(ds, r_mc, 1);

    SharedDataStream sds(ds);

    sds << mc.ensmble << mc.rangenerator
        << mc.naccept << mc.nreject
        << static_cast<const Move&>(mc);

    return ds;
}

/** Deserialise from a binary data stream */
QDataStream SIREMOVE_EXPORT &operator>>(QDataStream &ds, MonteCarlo &mc)
{
    VersionID v = readHeader(ds, r_mc);

    if (v == 1)
    {
        SharedDataStream sds(ds);
    
        sds >> mc.ensmble
            >> mc.rangenerator
            >> mc.naccept >> mc.nreject
            >> static_cast<Move&>(mc);
    }
    else
        throw version_error(v, "1", r_mc, CODELOC);

    return ds;
}

/** Construct using the supplied random number generator */
MonteCarlo::MonteCarlo() : Move(), naccept(0), nreject(0)
{}

/** Copy constructor */
MonteCarlo::MonteCarlo(const MonteCarlo &other)
           : Move(other), ensmble(other.ensmble),
             rangenerator(other.rangenerator),
             naccept(other.naccept), nreject(other.nreject)
{}

/** Destructor */
MonteCarlo::~MonteCarlo()
{}

/** Internal function called by derived classes to set the ensemble
    for this move */
void MonteCarlo::setEnsemble(const Ensemble &ensemble)
{
    ensmble = ensemble;
}

/** Return the ensemble for this move */
Ensemble MonteCarlo::ensemble() const
{
    return ensmble;
}

/** Copy assignment */
MonteCarlo& MonteCarlo::operator=(const MonteCarlo &other)
{
    ensmble = other.ensmble;
    rangenerator = other.rangenerator;
    naccept = other.naccept;
    nreject = other.nreject;

    Move::operator=(other);

    return *this;
}

/** Comparison operator */
bool MonteCarlo::operator==(const MonteCarlo &other) const
{
    return rangenerator == other.rangenerator and
           naccept == other.naccept and
           nreject == other.nreject;
}

/** Comparison operator */
bool MonteCarlo::operator!=(const MonteCarlo &other) const
{
    return not this->operator==(other);
}

/** Set the random number generator to use for these moves */
void MonteCarlo::setGenerator(const RanGenerator &generator)
{
    rangenerator = generator;
}

/** Return the random number generator used for these moves */
const RanGenerator& MonteCarlo::generator() const
{
    return rangenerator;
}

/** Return the number of attempted moves */
quint32 MonteCarlo::nAttempted() const
{
    return naccept + nreject;
}

/** Return the number of accepted moves */
quint32 MonteCarlo::nAccepted() const
{
    return naccept;
}

/** Return the number of rejected moves */
quint32 MonteCarlo::nRejected() const
{
    return nreject;
}

/** Return the acceptance ratio (ratio of
    accepted moves to attempted moves) */
double MonteCarlo::acceptanceRatio() const
{
    int ntotal = this->nAttempted();
    
    if (ntotal == 0)
        return 0;
    else
        return double(this->nAccepted()) / double(ntotal);
}

/** Zero the move statistics */
void MonteCarlo::clearStatistics()
{
    naccept = 0;
    nreject = 0;
}

/** Perform the NVT Monte Carlo test, using the supplied change in energy
    and the supplied change in biasing probabilities */
bool MonteCarlo::test(double new_energy, double old_energy,
                      double new_bias, double old_bias)
{
    double beta = -1.0 / (k_boltz * ensmble.temperature().value());

    double x = (new_bias / old_bias) * std::exp( beta*(new_energy - old_energy) );

    if (x > 1 or x > rangenerator.rand())
    {
        ++naccept;
        return true;
    }
    else
    {
        ++nreject;
        return false;
    }
}

/** Perform the NVT Monte Carlo test, using the supplied change in energy
    (no change in biasing factor) */
bool MonteCarlo::test(double new_energy, double old_energy)
{
    if (new_energy <= old_energy)
    {
        ++naccept;
        return true;
    }

    double beta = -1.0 / (k_boltz * ensmble.temperature().value());

    double x = std::exp( beta*(new_energy - old_energy) );

    if (x > 1 or x > rangenerator.rand())
    {
        ++naccept;
        return true;
    }
    else
    {
        ++nreject;
        return false;
    }
}

/** Perform the NPT Monte Carlo test, using the supplied change in 
    energy and supplied change in volume (no change in biasing factor) */
bool MonteCarlo::test(double new_energy, double old_energy,
                      int nmolecules,
                      const Volume &new_volume, const Volume &old_volume)
{
    double p_deltav = this->pressure() * (new_volume - old_volume);

    qDebug() << this->pressure() << (new_volume - old_volume);

    double vratio = nmolecules * std::log(new_volume / old_volume);

    double beta = -1.0 / (k_boltz * ensmble.temperature().value());

    qDebug() << "P dV =" << p_deltav << "vratio = " << vratio
             << "beta =" << beta << "dE = " << new_energy - old_energy;

    double x = std::exp( beta * (new_energy - old_energy + p_deltav) + vratio );

    qDebug() << "x =" << x << "delta =" << (-beta * (new_energy - old_energy + p_deltav) + vratio);

    if ( x > 1 or x > rangenerator.rand() )
    {
        ++naccept;
        return true;
    }
    else
    {
        ++nreject;
        return false;
    }
}

/** Perform the NPT Monte Carlo test, using the supplied change in 
    energy and supplied change in volume (with a change in biasing factor) */
bool MonteCarlo::test(double new_energy, double old_energy,
                      int nmolecules,
                      const Volume &new_volume, const Volume &old_volume,
                      double new_bias, double old_bias)
{
    double p_deltav = this->pressure() * (new_volume - old_volume);

    double vratio = nmolecules * std::log(new_volume / old_volume);

    double beta = -1.0 / (k_boltz * ensmble.temperature().value());

    double x =  (new_bias / old_bias) * 
                    std::exp( -beta * (new_energy - old_energy + p_deltav) + vratio );

    if ( x > 1 or x > rangenerator.rand() )
    {
        ++naccept;
        return true;
    }
    else
    {
        ++nreject;
        return false;
    }
}
