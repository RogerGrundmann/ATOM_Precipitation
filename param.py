# Given a parameter definition, generates necessary C++ bindings:
#   cube/CubeParameters.h   — member-variable declarations (included in cCubeModel class body)
#   cube/CubeParameters.cpp — SetDefaultConfig() implementation
# coding=utf-8

import os

def main():

    # Tuples: (name, description, datatype, default)
    # L_atm is already a static const in cCubeModel.h — omit it here.







    PARAMS = {
        'common': [
            ('output_path',             'directory where model outputs should be placed (must end in /)', 'string',  'output_ATOM_Cube_Turbulence/'),
            ('verbose',                 'enable verbose diagnostic output',                               'bool',    False),
#            ('paraview_panorama_vts_flag', 'flag to control if create paraview panorama',                 'bool',    False),
            ('paraview_panorama_vts_flag', 'flag to control if create paraview panorama',                 'bool',    True),


#            ('L_atm',                    'extension of the atmosphere shell in m, Benocci',                'double', 500.0),
            ('L_atm',                   'extension of the atmosphere shell in m, Benocci',                'double', 1.0),
#            ('L_atm',                   'extension of the atmosphere shell in m, Benocci',                'double', 0.2),

            ('p_0',                     'pressure at sea level in hPa',                                   'double',  1013.25),

            ('sigma',                   'Stefan-Boltzmann constant W/(m²*K4)',                            'double',  5.670280e-8),
            ('eps_residuum',            'relative error, end of iterations reached, 1% error allowed',    'double',  1.0e-4),

#            ('turb_model',             'turbulence model: none, k_epsilon, k_omega, k_omega_SST',         'string',  'none'),
            ('turb_model',              'turbulence model: none, k_epsilon, k_omega, k_omega_SST',        'string',  'k_omega_SST'),
#            ('turb_model',             'turbulence model: none, k_epsilon, k_omega, k_omega_SST',         'string',  'k_omega'),
#            ('turb_model',             'turbulence model: none, k_epsilon, k_omega, k_omega_SST',         'string',  'k_epsilon'),
        ],

        'cube': [
#            ('nm',            'the maximum number of iterations',                                         'int',     4),
#            ('nm',            'the maximum number of iterations',                                         'int',     60),
#            ('nm',            'the maximum number of iterations',                                         'int',     120),
#            ('nm',            'the maximum number of iterations',                                         'int',     256),
#            ('nm',            'the maximum number of iterations',                                         'int',     560),
            ('nm',            'the maximum number of iterations',                                         'int',     1120),
#            ('checkpoint',    'control when to write output files',                                       'int',     20),
#            ('checkpoint',    'control when to write output files',                                       'int',     56),
            ('checkpoint',    'control when to write output files',                                       'int',     112),
#            ('panorama_print','control when to write panorama files',                                     'int',     10),
            ('panorama_print','control when to write panorama files',                                     'int',     560),

            ('re',            'Reynolds number for laminar flows',                                        'double',  1000.0),
            ('pr',            'Prandtl number of air for laminar flows',                                  'double',  0.7179),
            ('pr_turb',       'turbulent Prandtl number for temperature transport',                       'double',  0.9),
            ('R_Air',         'specific gas constant of air in J/(kg*K)',                                 'double',  286.9),
#            ('u_0',           'annual mean of surface wind velocity in m/s',                              'double',  8.0),
            ('u_0',           'annual mean of surface wind velocity in m/s, Benocci',                     'double',  1.0),
#            ('u_0',           'annual mean of surface wind velocity in m/s',                              'double',  0.02),
#            ('u_0', 'annual mean of surface wind velocity in m/s, Benocci',                               'double', 5.81), 
    ]
}
    # Map param.py type names to C++ types
    cpp_type = {
        'string': 'std::string',
        'bool':   'bool',
        'int':    'int',
        'double': 'double',
        'float':  'float',
    }

    # Helper: format a default value as a C++ literal
    def cpp_default(dtype, val):
        if dtype == 'string':
            return f'"{val}"'
        if dtype == 'bool':
            return 'true' if val else 'false'
        if dtype == 'double':
            # Preserve scientific notation when appropriate
            return repr(float(val)).replace('e-0', 'e-').replace('e+0', 'e+')
        return str(val)

    # ------------------------------------------------------------------ #
    # Generate cube/CubeParameters.h
    # This file is #included inside the cCubeModel class body (public:).
    # It only contains member-variable declarations.
    # ------------------------------------------------------------------ #
    header_path = os.path.join('cube', 'CubeParameters.h')
    with open(header_path, 'w') as f:
        f.write('// cube/CubeParameters.h\n')
        f.write('// THIS FILE IS AUTOMATICALLY GENERATED BY param.py\n')
        f.write('// ANY CHANGES WILL BE OVERWRITTEN AT COMPILE TIME\n')
        f.write('\n')
        for section, params in PARAMS.items():
            f.write(f'    // {section} section\n')
            for name, desc, dtype, default in params:
                ctype = cpp_type[dtype]
                f.write(f'    {ctype} {name};  // {desc}\n')
            f.write('\n')

    # ------------------------------------------------------------------ #
    # Generate cube/CubeParameters.cpp
    # Implements cCubeModel::SetDefaultConfig().
    # ------------------------------------------------------------------ #
    cpp_path = os.path.join('cube', 'CubeParameters.cpp')
    with open(cpp_path, 'w') as f:
        f.write('// cube/CubeParameters.cpp\n')
        f.write('// THIS FILE IS AUTOMATICALLY GENERATED BY param.py\n')
        f.write('// ANY CHANGES WILL BE OVERWRITTEN AT COMPILE TIME\n')
        f.write('\n')
        f.write('#include "cCubeModel.h"\n')
        f.write('\n')
        f.write('void cCubeModel::SetDefaultConfig() {\n')
        f.write('\n')
        for section, params in PARAMS.items():
            f.write(f'    // {section} section\n')
            for name, desc, dtype, default in params:
                val = cpp_default(dtype, default)
                f.write(f'    {name} = {val};\n')
            f.write('\n')
        f.write('}\n')

    print(f'Generated {header_path}')
    print(f'Generated {cpp_path}')

    # ------------------------------------------------------------------ #
    # Generate python/cube.pxd
    # Cython extern declarations for cCubeModel, consumed by pycube.pyx.
    # ------------------------------------------------------------------ #
    cython_type = {
        'string': 'string',
        'bool':   'bool',
        'int':    'int',
        'double': 'double',
        'float':  'float',
    }

    pxd_path = os.path.join('python', 'cube.pxd')
    with open(pxd_path, 'w') as f:
        f.write('# python/cube.pxd\n')
        f.write('# THIS FILE IS AUTOMATICALLY GENERATED BY param.py\n')
        f.write('# ANY CHANGES WILL BE OVERWRITTEN AT COMPILE TIME\n')
        f.write('from libcpp cimport bool\n')
        f.write('from libcpp.string cimport string\n')
        f.write('from libcpp.vector cimport vector\n')
        f.write('\n')
        f.write('cdef extern from "cCubeModel.h":\n')
        f.write('    cppclass cCubeModel:\n')
        f.write('        cCubeModel() except +\n')
        f.write('        void LoadConfig(const char *filename)\n')
        f.write('        void Run()\n')
        f.write('        void RunTimeSlice(int time_slice)\n')
        f.write('        vector[float] get_layer_heights()\n')
        for section, params in PARAMS.items():
            f.write(f'        # {section} section\n')
            for name, desc, dtype, default in params:
                ctype = cython_type[dtype]
                f.write(f'        {ctype} {name}\n')

    print(f'Generated {pxd_path}')

    # ------------------------------------------------------------------ #
    # Generate python/pycube.pyx
    # Fill the {{ cube_params }} placeholder in the template with Cython
    # property wrappers for every parameter.
    # ------------------------------------------------------------------ #
    cython_cast = {
        'string': 'string',
        'bool':   'bool',
        'int':    'int',
        'double': 'double',
        'float':  'float',
    }

    props = []
    for section, params in PARAMS.items():
        props.append(f'    # {section} section')
        for name, desc, dtype, default in params:
            cast = cython_cast[dtype]
            props.append(f'    property {name}:')
            props.append(f'        def __get__(CubeModel self):')
            props.append(f'            self._check_alive()')
            props.append(f'            return self._thisptr.{name}')
            props.append(f'    ')
            props.append(f'        def __set__(CubeModel self, value):')
            props.append(f'            self._check_alive()')
            props.append(f'            self._thisptr.{name} = <{cast}> value')
            props.append(f'    ')
    cube_params_block = '\n'.join(props)

    template_path = os.path.join('python', 'pyatom.pyx.template')
    pyx_path      = os.path.join('python', 'pycube.pyx')
    with open(template_path) as f:
        template = f.read()
    pyx_content = template.replace('{{ cube_params }}', cube_params_block)
    # Prepend the standard auto-generated header comment
    header = (
        '# python/pycube.pyx\n'
        '# THIS FILE IS AUTOMATICALLY GENERATED BY param.py\n'
        '# ANY CHANGES WILL BE OVERWRITTEN AT COMPILE TIME\n'
        '\n'
    )
    with open(pyx_path, 'w') as f:
        f.write(header + pyx_content)

    print(f'Generated {pyx_path}')


if __name__ == '__main__':
    main()
