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

     # MPP unit tests — synthetic / generated data only. Tests that compared
     # against the PSS Python reference oracle (which lives outside this
     # repo) have been removed; everything still here either builds its own
     # fixture in-process or asserts pure-math invariants. No env vars
     # required.
     foreach t : ['mpp_rank_test', 'mpp_align_test', 'mpp_ap_test',
                  'mpp_shift_test', 'mpp_stack_test',
                  'mpp_pipeline_test', 'mpp_derot_geom_test']
         mpp_exec = executable(t,
                               t + '.cpp',
                               dependencies : [siril_dep, criterion_dep],
                               link_args : [siril_link_arg, '-Wl,--unresolved-symbols=ignore-all'],
                               c_args : siril_c_flag,
                               cpp_args : siril_cpp_flag)
         test(t, mpp_exec, suite: 'mpp')
     endforeach

     mpp_sidecar_exec = executable('mpp_sidecar_test',
                                   'mpp_sidecar_test.c',
                                   dependencies : [siril_dep, criterion_dep],
                                   link_args : [siril_link_arg, '-Wl,--unresolved-symbols=ignore-all'],
                                   c_args : siril_c_flag,
                                   cpp_args : siril_cpp_flag)

     test('mpp_sidecar_test', mpp_sidecar_exec, suite: 'mpp')

     mpp_derot_sidecar_exec = executable('mpp_derot_sidecar_test',
                                         'mpp_derot_sidecar_test.c',
                                         dependencies : [siril_dep, criterion_dep],
                                         link_args : [siril_link_arg, '-Wl,--unresolved-symbols=ignore-all'],
                                         c_args : siril_c_flag,
                                         cpp_args : siril_cpp_flag)

     test('mpp_derot_sidecar_test', mpp_derot_sidecar_exec, suite: 'mpp')

     planet_ephem_exec = executable('planet_ephem_test',
                                    'planet_ephem_test.c',
                                    dependencies : [siril_dep, criterion_dep],
                                    link_args : [siril_link_arg, '-Wl,--unresolved-symbols=ignore-all'],
                                    c_args : siril_c_flag,
                                    cpp_args : siril_cpp_flag)

     test('planet_ephem_test', planet_ephem_exec, suite: 'mpp')

     harmonize_img_t_exec = executable('harmonize_img_t_test',
                                       'harmonize_img_t_test.cpp',
                                       dependencies : [siril_dep, criterion_dep],
                                       link_args : [siril_link_arg, '-Wl,--unresolved-symbols=ignore-all'],
                                       c_args : siril_c_flag,
                                       cpp_args : siril_cpp_flag)

     test('harmonize_img_t_test', harmonize_img_t_exec, suite: 'harmonize')

     da3d_dft_perf_exec = executable('da3d_dft_perf',
                                     'da3d_dft_perf.cpp',
                                     dependencies : [siril_dep, criterion_dep],
                                     link_args : [siril_link_arg, '-Wl,--unresolved-symbols=ignore-all'],
                                     c_args : siril_c_flag,
                                     cpp_args : siril_cpp_flag)

     test('da3d_dft_perf', da3d_dft_perf_exec, suite: 'perfs', protocol: 'exitcode', is_parallel : false)

     # Performance testing
     sorting_perf_exec = executable('sorting_perf',
                                   'measure_sorting.c',
                                   dependencies : siril_dep,
                                   link_args : siril_link_arg,
                                   c_args : siril_c_flag,
                                   cpp_args : siril_cpp_flag)

     test('sorting_perf', sorting_perf_exec, suite: 'perfs', protocol: 'exitcode', is_parallel : false)

endif

