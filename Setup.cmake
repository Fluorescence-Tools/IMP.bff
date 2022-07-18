# libmongoc
###########################
FIND_PACKAGE (mongoc-1.0 1.7 REQUIRED)
LINK_LIBRARIES (mongo::mongoc_shared)
