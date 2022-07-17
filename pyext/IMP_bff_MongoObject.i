
// IMP_SWIG_BASE_OBJECT(IMP::bff, MongoObject, MongoObjects);

%template(ListMongoObjectPtr) std::list<IMP::bff::MongoObject*>;
%attributestring(IMP::bff::MongoObject, std::string, name, get_name, set_name);
%attributestring(IMP::bff::MongoObject, std::string, oid, get_own_oid, set_own_oid);
%attribute(IMP::bff::MongoObject, bool, is_connected_to_db, is_connected_to_db);

%shared_ptr(IMP::bff::MongoObject)

%include "IMP/bff/MongoObject.h"

%extend IMP::bff::MongoObject {

    public:
        %template(get_array_double) get_array<double>;
        %template(set_array_double) set_array<double>;
        %template(get_array_long) get_array<long>;
        %template(set_array_long) set_array<long>;
        %template(get_array_int) get_array<int>;
        %template(set_array_int) set_array<int>;

        %template(get_bool) get_<bool>;
        %template(set_bool) set_<bool>;
        %template(get_double) get_<double>;
        %template(set_double) set_<double>;
        %template(get_int) get_<int>;
        %template(set_int) set_<int>;

        %template(connect_object_to_db_mongo) connect_object_to_db<IMP::bff::MongoObject*>;

        IMP::bff::MongoObject* __getitem__(std::string key) {
            return(*self)[key];
        }

}

