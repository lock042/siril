#
# Siril Tests
#

# Criterion based unit-tests
if criterion_dep.found()
     sorting_exec = executable('sorting',
                               'sorting.c',
                               dependencies : [siril_dep, criterion_dep],
                               link_args : siril_link_arg,
                               c_args : siril_c_flag,
                               cpp_args : siril_cpp_flag)

     test('sorting', sorting_exec)

     stacking_exec = executable('stacking_blocks_test',
                                'stacking_blocks_test.c',
                                dependencies : [siril_dep, criterion_dep],
                                link_args : [siril_link_arg, '-Wl,--unresolved-symbols=ignore-all'],
                                c_args : siril_c_flag,
                                cpp_args : siril_cpp_flag)

     test('stacking_blocks_test', stacking_exec)
     
     date_exec = executable('siril_date_test',
                            'siril_date_test.c',
                            dependencies : [siril_dep, criterion_dep],
                            link_args : [siril_link_arg, '-Wl,--unresolved-symbols=ignore-all'],
                            c_args : siril_c_flag,
                            cpp_args : siril_cpp_flag)

     test('siril_date_test', date_exec)
     
     world_cs_exec = executable('siril_world_cs_test',
                                'siril_world_cs_test.c',
                                dependencies : [siril_dep, criterion_dep],
                                link_args : [siril_link_arg, '-Wl,--unresolved-symbols=ignore-all'],
                                c_args : siril_c_flag,
                                cpp_args : siril_cpp_flag)

     test('siril_world_cs_test', world_cs_exec)


     imoper_exec = executable('imoper_test',
                              'imoper_test.c',
                              dependencies : [siril_dep, criterion_dep],
                              link_args : [siril_link_arg, '-Wl,--unresolved-symbols=ignore-all'],
                              c_args : siril_c_flag,
                              cpp_args : siril_cpp_flag)

     test('imoper_test', imoper_exec, suite: 'arithmetic')

     soper_exec = executable('soper_test',
                             'soper_test.c',
                             dependencies : [siril_dep, criterion_dep],
                             link_args : [siril_link_arg, '-Wl,--unresolved-symbols=ignore-all'],
                             c_args : siril_c_flag,
                             cpp_args : siril_cpp_flag)

     test('soper_test', soper_exec, suite: 'arithmetic')

     photometry_exec = executable('photometry_test',
                                  'photometry_test.c',
                                  dependencies : [siril_dep, criterion_dep],
                                  link_args : [siril_link_arg, '-Wl,--unresolved-symbols=ignore-all'],
                                  c_args : siril_c_flag,
                                  cpp_args : siril_cpp_flag)

     test('photometry_test', photometry_exec, suite: 'science')

     rejection_exec = executable('rejection_test',
                                 'rejection_test.c',
                                 dependencies : [siril_dep, criterion_dep],
                                 link_args : [siril_link_arg, '-Wl,--unresolved-symbols=ignore-all'],
                                 c_args : siril_c_flag,
                                 cpp_args : siril_cpp_flag)

     test('rejection_test', rejection_exec, suite: 'arithmetic')

     ser_exec = executable('ser_test',
                           'ser_test.c',
                           dependencies : [siril_dep, criterion_dep],
                           link_args : [siril_link_arg, '-Wl,--unresolved-symbols=ignore-all'],
                           c_args : siril_c_flag,
                           cpp_args : siril_cpp_flag)

     test('ser_test', ser_exec)

     siril_spawn_test_exec = executable('siril_spawn_test',
                                        'siril_spawn_test.c',
                                        dependencies : [siril_dep, criterion_dep],
                                        link_args : siril_link_arg,
                                        c_args : siril_c_flag,
                                        cpp_args : siril_cpp_flag)

     test('siril_spawn_test', siril_spawn_test_exec)

     # ----- FLIS tests (stage 1.7) -----
     test_flis_lifecycle_exec = executable('test_flis_lifecycle',
                                'test_flis_lifecycle.c',
                                dependencies : [siril_dep, criterion_dep],
                                link_args : [siril_link_arg, '-Wl,--unresolved-symbols=ignore-all'],
                                c_args : siril_c_flag,
                                cpp_args : siril_cpp_flag)
     test('test_flis_lifecycle', test_flis_lifecycle_exec, suite: 'flis')

     test_flis_roundtrip_exec = executable('test_flis_roundtrip',
                                'test_flis_roundtrip.c',
                                dependencies : [siril_dep, criterion_dep],
                                link_args : [siril_link_arg, '-Wl,--unresolved-symbols=ignore-all'],
                                c_args : siril_c_flag,
                                cpp_args : siril_cpp_flag)
     test('test_flis_roundtrip', test_flis_roundtrip_exec, suite: 'flis')

     test_flis_geometry_exec = executable('test_flis_geometry',
                                'test_flis_geometry.c',
                                dependencies : [siril_dep, criterion_dep],
                                link_args : [siril_link_arg, '-Wl,--unresolved-symbols=ignore-all'],
                                c_args : siril_c_flag,
                                cpp_args : siril_cpp_flag)
     test('test_flis_geometry', test_flis_geometry_exec, suite: 'flis')

     test_flis_compose_exec = executable('test_flis_compose',
                                'test_flis_compose.c',
                                dependencies : [siril_dep, criterion_dep],
                                link_args : [siril_link_arg, '-Wl,--unresolved-symbols=ignore-all'],
                                c_args : siril_c_flag,
                                cpp_args : siril_cpp_flag)
     test('test_flis_compose', test_flis_compose_exec, suite: 'flis')

     test_flis_flatten_exec = executable('test_flis_flatten',
                                'test_flis_flatten.c',
                                dependencies : [siril_dep, criterion_dep],
                                link_args : [siril_link_arg, '-Wl,--unresolved-symbols=ignore-all'],
                                c_args : siril_c_flag,
                                cpp_args : siril_cpp_flag)
     test('test_flis_flatten', test_flis_flatten_exec, suite: 'flis')

     test_flis_undo_exec = executable('test_flis_undo',
                                'test_flis_undo.c',
                                dependencies : [siril_dep, criterion_dep],
                                link_args : [siril_link_arg, '-Wl,--unresolved-symbols=ignore-all'],
                                c_args : siril_c_flag,
                                cpp_args : siril_cpp_flag)
     test('test_flis_undo', test_flis_undo_exec, suite: 'flis')

     test_flis_worker_exec = executable('test_flis_worker',
                                'test_flis_worker.c',
                                dependencies : [siril_dep, criterion_dep],
                                link_args : [siril_link_arg, '-Wl,--unresolved-symbols=ignore-all'],
                                c_args : siril_c_flag,
                                cpp_args : siril_cpp_flag)
     test('test_flis_worker', test_flis_worker_exec, suite: 'flis')

     test_flis_cmd_exec = executable('test_flis_cmd',
                                'test_flis_cmd.c',
                                dependencies : [siril_dep, criterion_dep],
                                link_args : [siril_link_arg, '-Wl,--unresolved-symbols=ignore-all'],
                                c_args : siril_c_flag,
                                cpp_args : siril_cpp_flag)
     test('test_flis_cmd', test_flis_cmd_exec, suite: 'flis')

     # Performance testing
     sorting_perf_exec = executable('sorting_perf',
                                   'measure_sorting.c',
                                   dependencies : siril_dep,
                                   link_args : siril_link_arg,
                                   c_args : siril_c_flag,
                                   cpp_args : siril_cpp_flag)

     test('sorting_perf', sorting_perf_exec, suite: 'perfs', protocol: 'exitcode', is_parallel : false)

endif

