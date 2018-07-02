import fmodpy

test = fmodpy.fimport("test_alloc.f90", verbose=True, mod_name="test",
                      module_link_args=["-lblas","-llapack","-lgfortran"])

test.allocate_max_lapack_work(20, 40)
