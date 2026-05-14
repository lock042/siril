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

     # MPP / PSS port skeleton — see pss_port_plan.md. Real assertions land per phase.
     # MPP_PSS_TEST_DATA_DIR points at the bundled PSS reference oracle so
     # tests that compare against pre-computed PSS artifacts can find them;
     # if the dir is missing those tests skip themselves rather than failing.
     mpp_test_env = environment()
     mpp_test_env.set('MPP_PSS_TEST_DATA_DIR',
                      meson.project_source_root() / 'tools' / 'pss_reference')
     foreach t : ['mpp_rank_test', 'mpp_align_test', 'mpp_ap_test',
                  'mpp_shift_test', 'mpp_stack_test', 'mpp_drizzle_test',
                  'mpp_pipeline_test']
         mpp_exec = executable(t,
                               t + '.cpp',
                               dependencies : [siril_dep, criterion_dep],
                               link_args : [siril_link_arg, '-Wl,--unresolved-symbols=ignore-all'],
                               c_args : siril_c_flag,
                               cpp_args : siril_cpp_flag)
         test(t, mpp_exec, suite: 'mpp', env: mpp_test_env)
     endforeach

     mpp_sidecar_exec = executable('mpp_sidecar_test',
                                   'mpp_sidecar_test.c',
                                   dependencies : [siril_dep, criterion_dep],
                                   link_args : [siril_link_arg, '-Wl,--unresolved-symbols=ignore-all'],
                                   c_args : siril_c_flag,
                                   cpp_args : siril_cpp_flag)

     test('mpp_sidecar_test', mpp_sidecar_exec, suite: 'mpp', env: mpp_test_env)

     # Performance testing
     sorting_perf_exec = executable('sorting_perf',
                                   'measure_sorting.c',
                                   dependencies : siril_dep,
                                   link_args : siril_link_arg,
                                   c_args : siril_c_flag,
                                   cpp_args : siril_cpp_flag)

     test('sorting_perf', sorting_perf_exec, suite: 'perfs', protocol: 'exitcode', is_parallel : false)

endif

