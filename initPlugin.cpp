#include <string.h>
#include <stdio.h>
#include <sofa/helper/system/FileRepository.h>
#include <sofa/helper/system/SetDirectory.h>
#include <sofa/core/ObjectFactory.h>

#define Q(x) #x
#define QUOTE(x) Q(x)

namespace sofa {

namespace msofaplugin {

	//Here are just several convenient functions to help user to know what contains the plugin

	extern "C" {
        void initExternalModule();
        const char* getModuleName();
        const char* getModuleVersion();
        const char* getModuleLicense();
        const char* getModuleDescription();
        const char* getModuleComponentList();
	}
	
	void initExternalModule()
	{
		static bool first = true;
		if (first)
		{
            first = false;
#ifdef PLUGIN_DATA_DIR
            sofa::helper::system::DataRepository.addLastPath(std::string(QUOTE(PLUGIN_DATA_DIR)));
#endif
            sofa::helper::system::DataRepository.addLastPath(sofa::helper::system::SetDirectory::GetCurrentDir());
		}
	}

	const char* getModuleName()
	{
        return "SofaDefaultPlugin";
	}

	const char* getModuleVersion()
	{
        return "0";
	}

	const char* getModuleLicense()
	{
		return "LGPL";
	}

	const char* getModuleDescription()
    {
        return "Solver plugin";
	}

	const char* getModuleComponentList()
	{
        return "";
	}

} 

} 
