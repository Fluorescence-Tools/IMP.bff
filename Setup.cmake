# libmongoc
###########################
FIND_PACKAGE (mongoc-1.0 1.7 REQUIRED)
LINK_LIBRARIES (mongo::mongoc_shared)

# RTTR
###########################
FIND_PACKAGE(RTTR CONFIG REQUIRED Core)
LINK_LIBRARIES(RTTR::Core)
