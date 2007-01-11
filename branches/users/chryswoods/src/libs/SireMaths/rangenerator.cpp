
#include <QMutex>
#include <QVector>

#include <boost/noncopyable.hpp>
#include <boost/scoped_array.hpp>

#include <limits>

#include "rangenerator.h"
#include "vector.h"

#include "ThirdParty/MersenneTwister.h"

#include "SireError/errors.h"
#include "SireStream/datastream.h"
#include "SireStream/shareddatastream.h"

using namespace SireMaths;
using namespace SireMaths::detail;

using namespace SireStream;

namespace SireMaths
{
namespace detail
{
class RanGeneratorPvt;
}
}

QDataStream& operator<<(QDataStream&, const SireMaths::detail::RanGeneratorPvt&);
QDataStream& operator>>(QDataStream&, SireMaths::detail::RanGeneratorPvt&);

namespace SireMaths
{

namespace detail
{

/** This class provides the private implementation of RanGenerator.
    This class is explicitly shared by instances of RanGenerator,
    and is therefore thread-safe.
    
    @author Christopher Woods
*/
class RanGeneratorPvt : public boost::noncopyable
{

friend QDataStream& ::operator<<(QDataStream&, const RanGeneratorPvt&);
friend QDataStream& ::operator>>(QDataStream&, RanGeneratorPvt&);

typedef MTRand::uint32 MTUInt32;

public:
    /** Construct a generator with a random seed */
    RanGeneratorPvt() : boost::noncopyable(),
                        mutex(QMutex::NonRecursive), 
                        mersenne_generator()
    {}
    
    /** Construct a generator with a specified seed */
    RanGeneratorPvt(quint32 seed) : boost::noncopyable(),
                                    mutex(QMutex::NonRecursive),
                                    mersenne_generator(seed)
    {}
    
    /** Construct a generator with a specified seed */
    RanGeneratorPvt(const QVector<quint32> &s) : boost::noncopyable(), 
                                                 mutex(QMutex::NonRecursive),
                                                 mersenne_generator()
    {
        this->seed(s);
    }

    /** Destructor */
    ~RanGeneratorPvt()
    {}
    
    /** Return a random double in [0,1) */
    double rand()
    {
        QMutexLocker lkr(&mutex);
        return mersenne_generator.rand();
    }
    
    /** Fill the passed vector with random numbers */
    void rand(QVector<double> &array)
    {
        if (not array.isEmpty())
        {
            QMutexLocker lkr(&mutex);
        
            int sz = array.count();
            double *a = array.data();
        
            for (int i=0; i<sz; ++i)
                a[i] = mersenne_generator.rand();
        }
    }
    
    /** Return a random double in [0,1) which has a full
        53 bits of randomness (as opposed to 32bits) */
    double rand53()
    {
        QMutexLocker lkr(&mutex);
        return mersenne_generator.rand53();
    }
    
    /** Fill the passed vector with random numbers */
    void rand53(QVector<double> &array)
    {
        if (not array.isEmpty())
        {
            QMutexLocker lkr(&mutex);
        
            int sz = array.count();
            double *a = array.data();
        
            for (int i=0; i<sz; ++i)
                a[i] = mersenne_generator.rand53();
        }
    }
    
private:
    /** Return a random number from a normal (gaussian) 
        distribution with given mean and variance */
    double _randNorm(double mean, double variance)
    {
        double rand0 = mersenne_generator.rand53();
        double rand1 = mersenne_generator.rand53();
        
        // Return a real number from a normal (Gaussian) distribution with given
        // mean and variance by Box-Muller method
        double r = std::sqrt( -2.0 * log(1.0-rand0) ) * variance;
        double phi = 2.0 * 3.14159265358979323846264338328 * rand1;
        
        return mean + r * std::cos(phi);
    }
public:
    /** Return a random number from a normal (gaussian) 
        distribution with given mean and variance */
    double randNorm(double mean, double variance)
    {
        QMutexLocker lkr(&mutex);
        
        return _randNorm(mean,variance);
    }
    
    /** Fill the vector with random normal numbers with
        given mean and variance */
    void randNorm(QVector<double> &array, double mean, double variance)
    {
        if (not array.isEmpty())
        {
            QMutexLocker lkr(&mutex);
        
            int sz = array.count();
            double *a = array.data();
        
            for (int i=0; i<sz; ++i)
                a[i] = _randNorm(mean,variance);
        }
    }
        
private:    
    /** Return a random vector on a sphere of radius 1.0 */
    Vector _vectorOnSphere()
    {
        while( true )
        {
            //use von Neumann acceptance/rejection method
            Vector v;
        
            v.setX( 1.0 - 2.0 * mersenne_generator.rand53() );
            v.setY( 1.0 - 2.0 * mersenne_generator.rand53() );
            v.setZ( 1.0 - 2.0 * mersenne_generator.rand53() );
        
            double lgth2 = v.length2();
        
            if (lgth2 < 1)
            {
                v /= std::sqrt(lgth2);
                return v;
            }
        }
    }
    
public:
    /** Return a random vector on a sphere of radius 1.0 */
    Vector vectorOnSphere()
    {
        QMutexLocker lkr(&mutex);
        
        return _vectorOnSphere();
    }
    
    /** Fill the array with random vectors on a sphere of radius 1.0 */
    void vectorOnSphere(QVector<Vector> &array)
    {
        if (not array.isEmpty())
        {
            QMutexLocker lkr(&mutex);
        
            int sz = array.count();
            Vector *a = array.data();
        
            for (int i=0; i<sz; ++i)
                a[i] = _vectorOnSphere();
        }
    }
    
    /** Return random 32bit unsigned integer (0 -> 2^32-1) */
    quint32 randInt()
    {
        QMutexLocker lkr(&mutex);
        return mersenne_generator.randInt();
    }
    
    /** Fill the array with random 32bit unsigned integers */
    void randInt(QVector<quint32> &array)
    {
        if (not array.isEmpty())
        {
            QMutexLocker lkr(&mutex);
        
            int sz = array.count();
            quint32 *a = array.data();
        
            for (int i=0; i<sz; ++i)
                a[i] = mersenne_generator.randInt();
        }
    }
    
    /** Return random 32bit unsigned integer (0 -> maxval) */
    quint32 randInt(quint32 maxval)
    {
        QMutexLocker lkr(&mutex);
        return mersenne_generator.randInt(maxval);
    }
    
    /** Fill the array with random 32bit unsigned integers */
    void randInt(QVector<quint32> &array, quint32 maxval)
    {
        if (not array.isEmpty())
        {
            QMutexLocker lkr(&mutex);
        
            int sz = array.count();
            quint32 *a = array.data();
        
            MTRand::uint32 m = maxval;
        
            for (int i=0; i<sz; ++i)
                a[i] = mersenne_generator.randInt(m);
        }
    }
    
private:
    /** Return random 64bit unsigned integer (0 -> 2^64-1) */
    quint64 _randInt64()
    {
        quint64 ran0 = mersenne_generator.randInt();
        quint64 ran1 = mersenne_generator.randInt();
        
        return (ran0 << 32) | ran1;
    }

    /** Return random 64bit unsigned integer (0 -> maxval) */
    quint64 _randInt64(quint64 maxval)
    {
        if (maxval <= std::numeric_limits<quint32>::max())
            //maxval can fit into a 32bit int - there is no
            //point using a 64bit generator!
            return mersenne_generator.randInt( quint32(maxval) );
    
        //use same algorithm in MersenneTwister.h
        quint64 used = maxval;
        
        used |= used >> 1;
        used |= used >> 2;
        used |= used >> 4;
        used |= used >> 8;
        used |= used >> 16;
        used |= used >> 32;
        
        quint64 i;
        
        do
        {
            i = _randInt64() & used;
        } while ( i > maxval );
        
        return i;
    }

public:
    /** Return random 64bit unsigned integer (0 -> 2^64-1) */
    quint64 randInt64()
    {
        QMutexLocker lkr(&mutex);
        return _randInt64();
    }
    
    /** Fill the array with random 64bit unsigned integers */
    void randInt64(QVector<quint64> &array)
    {
        if (not array.isEmpty())
        {
            QMutexLocker lkr(&mutex);
        
            int sz = array.count();
            quint64 *a = array.data();
        
            for (int i=0; i<sz; ++i)
                a[i] = _randInt64();
        }
    }
    
    /** Return random 64bit unsigned integer (0 -> maxval) */
    quint64 randInt64(quint64 maxval)
    {
        QMutexLocker lkr(&mutex);

        return _randInt64(maxval);
    }
    
    /** Fill the array with random 64bit unsigned integers */
    void randInt64(QVector<quint64> &array, quint64 maxval)
    {
        if (not array.isEmpty())
        {
            QMutexLocker lkr(&mutex);
        
            int sz = array.count();
            quint64 *a = array.data();
        
            for (int i=0; i<sz; ++i)
                a[i] = _randInt64(maxval);
        }
    }
    
    /** Randomly seed the generator */
    void seed()
    {
        QMutexLocker lkr(&mutex);
        
        mersenne_generator.seed();
    }
    
    /** Reseed the generator from an array of uints */
    void seed(const QVector<quint32> &s)
    {
        if (s.isEmpty())
        {
            this->seed();
            return;
        }
    
        QMutexLocker lkr(&mutex);
    
        //we need to convert this into an array of MTRand::uint32
        int sz = s.count();
            
        boost::scoped_array<MTUInt32> array( new MTUInt32[sz] );
            
        for (int i=0; i<sz; ++i)
            array[i] = s.constData()[i];
                
        mersenne_generator.seed( array.get(),  sz );
    }
    
    /** Reseed the generator */
    void seed(quint32 s)
    {
        QMutexLocker lkr(&mutex);
        
        mersenne_generator.seed(s);
    }
    
    /** Return an array containing the state of the random generator */
    QVector<quint32> getState()
    {
        mutex.lock();
        
        //create an array to hold the state
        boost::scoped_array<MTUInt32> array( new MTUInt32[MTRand::N] );
        
        mersenne_generator.save(array.get());
        
        mutex.unlock();
        
        //copy the array into a QVector<quint32>
        QVector<quint32> ret(MTRand::N);
        
        quint32 *ret_array = ret.data();
        
        for (int i=0; i<MTRand::N; ++i)
            ret_array[i] = array[i];
            
        return ret;
    }
    
    /** Load the state from an array - the array must have size
        MTRand::N (624) 
        
        \throw SireError::incompatible_error
    */
    void loadState(const QVector<quint32> &state)
    {
        //check that the array is of the right size...
        if (state.count() != MTRand::N)
            throw SireError::incompatible_error( QObject::tr(
                "Can only restore the state from an array of size %1, "
                "while you have provided an array of size %2.")
                    .arg(MTRand::N).arg(state.count()), CODELOC );
    
        //convert the array to the right type...
        boost::scoped_array<MTUInt32> array( new MTUInt32[MTRand::N] );
        
        const quint32 *state_array = state.constData();
        
        for (int i=0; i<MTRand::N; ++i)
            array[i] = state_array[i];
            
        QMutexLocker lkr(&mutex);
        
        mersenne_generator.load(array.get());
    }
    
private:
    /** Mutex to serialise access to the generator */
    QMutex mutex;

    /** The actual generator (Mersenne Twister) */
    MTRand mersenne_generator;
};

} // end of namespace detail

} // end of namespace SireMaths

////////////
//////////// Implementation of RanGeneratorPvt
////////////

static const RegisterMetaType<RanGeneratorPvt> r_ranpvt(MAGIC_ONLY,
                                  "SireMaths::detail::RanGeneratorPvt");

/** Serialise to a binary datastream */
QDataStream SIREMATHS_EXPORT &operator<<(QDataStream &ds, const RanGeneratorPvt &ranpvt)
{
    writeHeader(ds, r_ranpvt, 1);

    //this is a large object - full state is 624 quint32 (~2.5kB)
    typedef MTRand::uint32 MTUInt32;
    
    boost::scoped_array<MTUInt32> array( new MTUInt32[MTRand::N] );
    
    const_cast<QMutex&>(ranpvt.mutex).lock();
    ranpvt.mersenne_generator.save(array.get());
    const_cast<QMutex&>(ranpvt.mutex).unlock();
    
    ds << quint32(MTRand::N);
    
    for (int i=0; i<MTRand::N; ++i)
        //convert to platform neutral quint32
        ds << quint32( array[i] );
    
    return ds;
}

/** Deserialise from a binary datastream */
QDataStream SIREMATHS_EXPORT &operator>>(QDataStream &ds, RanGeneratorPvt &ranpvt)
{
    VersionID v = readHeader(ds, r_ranpvt);

    if (v == 1)
    {
        quint32 n;
    
        ds >> n;
        
        throw SireError::program_bug( QObject::tr(
            "The size of the Mersenne Twister state array has changed without "
            "the version number changing - this is a bug! (%1 vs. %2)")
                .arg(n).arg(MTRand::N), CODELOC );
                
        //this is a large object - full state is 624 quint32 (~2.5kB)
        typedef MTRand::uint32 MTUInt32;
    
        boost::scoped_array<MTUInt32> array( new MTUInt32[MTRand::N] );
        
        for (int i=0; i<MTRand::N; ++i)
        {
            quint32 val;
            ds >> val;
            
            array[i] = val;
        }
        
        QMutexLocker lkr( &(ranpvt.mutex) );
        ranpvt.mersenne_generator.load(array.get());
        ranpvt.mutex.unlock();
    }
    else
        throw version_error(v, "1", r_ranpvt, CODELOC);

    return ds;
}

////////////
//////////// Implementation of RanGenerator
////////////

static const RegisterMetaType<RanGenerator> r_rangen;

/** Serialise to a binary data stream */
QDataStream SIREMATHS_EXPORT &operator<<(QDataStream &ds, const RanGenerator &rangen)
{
    writeHeader(ds, r_rangen, 1);
    
    //SharedDataStream sds(ds) << rangen.d;
    ds << *(rangen.d);
    
    return ds;
}

/** Deserialise from a binary data stream */
QDataStream SIREMATHS_EXPORT &operator>>(QDataStream &ds, RanGenerator &rangen)
{
    VersionID v = readHeader(ds, r_rangen);
    
    if (v == 1)
    {
        if (not rangen.d.unique())
            rangen.d.reset( new RanGeneratorPvt() );
    
        //SharedDataStream sds(ds) >> rangen.d;
        ds >> *(rangen.d);
    }
    else
        throw version_error(v, "1", r_rangen, CODELOC);
    
    return ds;
}

static boost::shared_ptr<RanGeneratorPvt> shared_null( new RanGeneratorPvt() );

/** Create a randomly seeded generator
    (actually a copy of the global, random generator) */
RanGenerator::RanGenerator() : d(shared_null)
{}

/** Create a generator seeded with 'seed' */
RanGenerator::RanGenerator(quint32 seed)
             : d( new RanGeneratorPvt(seed) )
{}

/** Create a generator seeded with 'seed' */
RanGenerator::RanGenerator(const QVector<quint32> &seed)
             : d( new RanGeneratorPvt(seed) )
{}

/** Copy constructor - this takes an explicitly shared
    copy of 'other' (this is to prevent repeat random numbers
    from being generated by implicit copies!) */
RanGenerator::RanGenerator(const RanGenerator &other)
             : d( other.d )
{}

/** Destructor */
RanGenerator::~RanGenerator()
{}

/** Copy assignment */
RanGenerator& RanGenerator::operator=(const RanGenerator &other)
{
    d = other.d;
    return *this;
}

/** See the generator with a new, random seed - this will detach
    this explicitly shared copy of the generator */
void RanGenerator::seed()
{
    if (d.unique())
        d->seed();
    else
        d.reset( new RanGeneratorPvt() );
}

/** Seed the generator with 's'  - this will detach
    this explicitly shared copy of the generator */
void RanGenerator::seed(quint32 s)
{
    if (d.unique())
        d->seed(s);
    else
        d.reset( new RanGeneratorPvt(s) );
}

/** Seed the generator with 'seed' - this will detach
    this explicitly shared copy of the generator */
void RanGenerator::seed(const QVector<quint32> &s)
{
    if (d.unique())
        d->seed(s);
    else
        d.reset( new RanGeneratorPvt(s) );
}

/** Seed the generator with another generator - this 
    really just copies the generator as they are 
    all explicit copies of one another! */
void RanGenerator::seed(const RanGenerator &other)
{
    d = other.d;
}

/** Return a random real number on [0,1) */
double RanGenerator::rand()
{
    return d->rand();
}

/** Return an array of random real numbers on [0,1).
    If n == 0 then fill the array, otherwise generate
    n numbers. */
void RanGenerator::rand(QVector<double> &array, uint n)
{
    if (n)
        array.resize(n);

    d->rand(array);
}

/** Return a high-precision random real number on [0,1) */
double RanGenerator::rand53()
{
    return d->rand53();
}

/** Return an array of high-precision random real numbers on [0,1)
    If n == 0 then fill the array, otherwise generate
    n numbers. */
void RanGenerator::rand53(QVector<double> &array, uint n)
{
    if (n)
        array.resize(n);
        
    d->rand53(array);
}

/** Return a high-precision random number from the normal distribution 
    with supplied mean and variance. */
double RanGenerator::randNorm(double mean, double variance)
{
    return d->randNorm(mean, variance);
}

/** Return an array of high-precision random numbers
    from the normal distribution with supplied mean and variance.
    If n == 0 then fill the array, otherwise generate
    n numbers. */
void RanGenerator::randNorm(QVector<double> &array, double mean,
                            double variance, uint n)
{
    if (n)
        array.resize(n);
        
    d->randNorm(array, mean, variance);
}

/** Return a random vector on the unit sphere */
Vector RanGenerator::vectorOnSphere()
{
    return d->vectorOnSphere();
}

/** Return an array of random vectors on the unit sphere.
    If n == 0 then fill the array, otherwise generate
    n vectors. */
void RanGenerator::vectorOnSphere(QVector<Vector> &array, uint n)
{
    if (n)
        array.resize(n);
        
    d->vectorOnSphere(array);
}

/** Return a random 32bit unsigned integer (0 -> 2^32 - 1) */
quint32 RanGenerator::randInt()
{
    return d->randInt();
}

/** Return a random 32bit unsigned integer (0 -> maxval) */
quint32 RanGenerator::randInt(quint32 maxval)
{
    return d->randInt(maxval);
}

/** Return an array of random 32bit unsigned integers.
    If n == 0 then fill the array, otherwise generate
    n numbers. */
void RanGenerator::randInt(QVector<quint32> &array, uint n)
{
    if (n)
        array.resize(n);
        
    d->randInt(array);
}

/** Return an array of random 32bit unsigned integers.
    If n == 0 then fill the array, otherwise generate
    n numbers. */
void RanGenerator::randInt(QVector<quint32> &array, quint32 maxval, uint n)
{
    if (n)
        array.resize(n);
        
    d->randInt(array, maxval);
}

/** Return a random 64bit unsigned integer (0 -> 2^64 - 1) */
quint64 RanGenerator::randInt64()
{
    return d->randInt64();
}

/** Return a random 64bit unsigned integer (0 -> maxval) */
quint64 RanGenerator::randInt64(quint64 maxval)
{
    return d->randInt64(maxval);
}

/** Return an array of random 64bit unsigned integers.
    If n == 0 then fill the array, otherwise generate
    n numbers. */
void RanGenerator::randInt64(QVector<quint64> &array, uint n)
{
    if (n)
        array.resize(n);
        
    d->randInt64(array);
}

/** Return an array of random 64bit unsigned integers.
    If n == 0 then fill the array, otherwise generate
    n numbers. */
void RanGenerator::randInt64(QVector<quint64> &array, quint64 maxval, uint n)
{
    if (n)
        array.resize(n);
        
    d->randInt64(array, maxval);
}

/** Return the current state of the random number generator.
    Use this if you truly wish to get reproducible sequences
    of random numbers */
QVector<quint32> RanGenerator::getState()
{
    return d->getState();
}

/** Load the state into this generator - the state must have
    been produced by the getState() function above.
    
    This will detach this copy from shared storage.
    
    \throw SireError::incompatible_error
*/
void RanGenerator::setState(const QVector<quint32> &state)
{
    if (d.unique())
    {
        d->loadState(state);
    }
    else
    {
        d.reset( new RanGeneratorPvt() );
        d->loadState(state);
    }
}
