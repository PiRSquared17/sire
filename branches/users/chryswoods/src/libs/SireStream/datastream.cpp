
#include <QDataStream>

#include "datastream.h"
#include "getmagic.h"
#include "magic_error.h"
#include "version_error.h"
#include "registermetatype.hpp"

namespace SireStream
{

/** Write a header to the data stream that describes the type and version
    of the object that is about to be written */
QDataStream SIRESTREAM_EXPORT &writeHeader(QDataStream &ds, 
                                           const RegisterMetaTypeBase &r_type,
                                           VersionID version)
{
    ds << r_type.magicID() << version;
    return ds;
}

/** Write a header to the data stream that contains the magic and version
    of the object that is about to be written */
QDataStream SIRESTREAM_EXPORT &writeHeader(QDataStream &ds, 
                                           MagicID magicid, 
                                           VersionID version)
{
    ds << magicid << version;
    return ds;
}

/** Read the header of the binary object to check that the type is correct
    and to obtain the binary data version */
VersionID SIRESTREAM_EXPORT readHeader(QDataStream &ds, 
                                       const RegisterMetaTypeBase &r_type)
{
    MagicID id;
    VersionID v;

    ds >> id >> v;

    if (id != r_type.magicID())
        throw SireStream::magic_error( id, r_type, CODELOC );

    return v;
}

/** Read the header of the binary object to check that the type is correct
    and to obtain the binary data version */
VersionID SIRESTREAM_EXPORT readHeader(QDataStream &ds, 
                                       MagicID magicid, 
                                       const char *type_name)
{
    MagicID id;
    VersionID v;

    ds >> id >> v;

    if (id != magicid)
        throw SireStream::magic_error( id, magicid, type_name, CODELOC );

    return v;
}

}
