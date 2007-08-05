
import os
import sys

from pyplusplus.module_builder import module_builder_t
from pyplusplus.module_builder import call_policies

from pygccxml.declarations.matchers import access_type_matcher_t
from pygccxml import declarations

sys.path.append("..")
from sireutils import *

wrap_classes = [ "AABox",
                 "Cartesian",
                 "CoordGroupBase",
                 "CoordGroup",
                 "CoordGroupEditor",
                 "MapFromCartesianFunction",
                 "MapFromSelfFunction",
                 "MappingFunctionBase",
                 "PeriodicBox",
                 "SpaceBase",
                 
                 "ConcreteProperty<SireVol::Cartesian, SireVol::SpaceBase>",
                 "ConcreteProperty<SireVol::PeriodicBox, SireVol::Cartesian>"
               ]

huge_classes = []

aliases = { "ConcreteProperty<SireVol::Cartesian, SireVol::SpaceBase>" :
            "ConcreteProperty_Cartesian_SpaceBase_",
            "ConcreteProperty<SireVol::PeriodicBox, SireVol::Cartesian>" :
            "ConcreteProperty_PeriodicBox_Cartesian_"}

extra_includes = [ "SireMaths/vector.h",
                   "SireMaths/quaternion.h",
                   "SireMaths/matrix.h" ]

def fix_coordgroupbase(c):
   c.decl("data").exclude()
   c.decl("constData").exclude()
   
   
def fix_coordgroup(c):   
   c.add_registration_code( "def( \"split\", &SireVol::CoordGroup::split< QVector<SireVol::CoordGroup> > )" )
   c.add_registration_code( "def( \"combine\", &SireVol::CoordGroup::combine< QVector<SireVol::CoordGroup> > )")
   
   c.add_registration_code("staticmethod(\"split\")")
   c.add_registration_code("staticmethod(\"combine\")")

def fix_coordgroupeditor(c):
   c.decl("data").exclude()
   
   c.decls( "translate" ).call_policies = call_policies.return_self()
   c.decls( "rotate" ).call_policies = call_policies.return_self()
   c.decls( "setCoordinates" ).call_policies = call_policies.return_self()
   

special_code = { "CoordGroupBase" : fix_coordgroupbase,
                 "CoordGroup" : fix_coordgroup,
                 "CoordGroupEditor" : fix_coordgroupeditor }

implicitly_convertible = [ ("QVector<SireMaths::Vector>","SireVol::CoordGroup"),
                           ("const SireVol::SpaceBase&","SireVol::Space"),
                           ("const SireVol::MappingFunctionBase&","SireVol::MappingFunction")
                         ]

incpaths = sys.argv[1:]
incpaths.insert(0, "../../")

#give the namespace and header files for the module to export
namespace = "SireVol"
modulename = "Vol"

headerfiles = ["sirevol_headers.h"]

#construct a module builder that will build the module's wrappers
mb = module_builder_t( files=headerfiles, 
                       include_paths=incpaths,
                       define_symbols=["SKIP_BROKEN_GCCXML_PARTS",
                                       "SKIP_TEMPLATE_DEFINITIONS"] )


populateNamespaces(mb)

for calldef in mb.calldefs():
    try:
      calldef.virtuality = declarations.VIRTUALITY_TYPES.NOT_VIRTUAL
    except:
      pass

#add calls to register hand-written wrappers
mb.add_declaration_code( "#include \"sirevol_containers.h\"" )
mb.add_registration_code( "register_SireVol_containers();", tail=False )

mb.calldefs().create_with_signature = True

#export each class in turn
for classname in wrap_classes:
   #tell the program to write wrappers for this class
   export_class(mb, classname, aliases, special_code)

register_implicit_conversions(mb, implicitly_convertible)

write_wrappers(mb, modulename, extra_includes, huge_classes)
