import IMP.bff as bff


DB_DICT = {
    'uri_string': "mongodb://localhost:27017",
    'db_string': "IMP.bff",
    'app_string': "bff",
    'collection_string': "test_collection"
}

def connects_to_db():
    mo = bff.MongoObject()
    mo.connect_to_db(**DB_DICT)
    return mo.is_connected_to_db 

CONNECTS = connects_to_db()
