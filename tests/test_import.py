import pytest

def test_namespace_alias():
    with pytest.raises(ImportError):
        from doranet import dn

def test_namespace_nesting():
    with pytest.raises(ImportError):
        from doranet import doranet