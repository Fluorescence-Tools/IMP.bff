"""The exception family is the same under both builds.

Under IMP, a C++ `IMP::ValueException` (or a `std::domain_error`, which IMP's
kernel maps to the same thing) reaches Python as `IMP.ValueException`, a
subclass of both `IMP.Exception` and `ValueError`. The standalone build has
no `IMP` package to put these on, so it makes the same classes itself
(`standalone/pyext/IMP_bff_standalone.macros.i`), and the IMP build exposes
them under `IMP.bff` as well (`pyext/swig.i-in`), so that a caller writes
`IMP.bff.ValueException` and is right in both worlds.

This test names nothing of IMP but `IMP.bff`, so it runs in both lanes.
"""

import unittest

import IMP.bff


class TestExceptionFamily(unittest.TestCase):

    def test_the_family_derives_as_imp_says(self):
        self.assertTrue(issubclass(IMP.bff.ValueException, IMP.bff.Exception))
        self.assertTrue(issubclass(IMP.bff.ValueException, ValueError))
        self.assertTrue(issubclass(IMP.bff.IOException, IOError))
        self.assertTrue(issubclass(IMP.bff.IndexException, IndexError))
        self.assertTrue(issubclass(IMP.bff.TypeException, TypeError))
        self.assertTrue(issubclass(IMP.bff.UsageException, IMP.bff.Exception))
        self.assertFalse(issubclass(IMP.bff.UsageException, ValueError))
        self.assertTrue(issubclass(IMP.bff.Exception, Exception))

    def test_a_refused_value_is_a_value_exception(self):
        # FitChiSquared::set_data throws std::domain_error; IMP maps that to
        # ValueException, and so does the standalone module.
        c = IMP.bff.FitChiSquared("chi2")
        with self.assertRaises(IMP.bff.ValueException):
            c.set_data([1.0, 2.0], [1.0])
        with self.assertRaises(ValueError):
            c.set_data([1.0, 2.0], [1.0])

    def test_a_missing_file_is_an_io_exception(self):
        with self.assertRaises(IMP.bff.IOException):
            IMP.bff.load_structure("/no/such/structure.pdb")
        with self.assertRaises(IOError):
            IMP.bff.load_structure("/no/such/structure.pdb")


if __name__ == "__main__":
    unittest.main()
