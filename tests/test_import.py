import pytest

def test_import_package():
    with pytest.raises(ImportError):
        import doranet as dn
