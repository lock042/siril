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
                                link_args : siril_link_arg,
                                c_args : siril_c_flag,
                                cpp_args : siril_cpp_flag)

     test('stacking_blocks_test', stacking_exec)

     date_exec = executable('siril_date_test',
                            'siril_date_test.c',
                            dependencies : [siril_dep, criterion_dep],
                            link_args : siril_link_arg,
                            c_args : siril_c_flag,
                            cpp_args : siril_cpp_flag)

     test('siril_date_test', date_exec)

     world_cs_exec = executable('siril_world_cs_test',
                                'siril_world_cs_test.c',
                                dependencies : [siril_dep, criterion_dep],
                                link_args : siril_link_arg,
                                c_args : siril_c_flag,
                                cpp_args : siril_cpp_flag)

     test('siril_world_cs_test', world_cs_exec)


     imoper_exec = executable('imoper_test',
                              'imoper_test.c',
                              dependencies : [siril_dep, criterion_dep],
                              link_args : siril_link_arg,
                              c_args : siril_c_flag,
                              cpp_args : siril_cpp_flag)

     test('imoper_test', imoper_exec, suite: 'arithmetic')

     fits_region_exec = executable('fits_region_test',
                                   'fits_region_test.c',
                                   dependencies : [siril_dep, criterion_dep],
                                   link_args : siril_link_arg,
                                   c_args : siril_c_flag,
                                   cpp_args : siril_cpp_flag)

     test('fits_region_test', fits_region_exec, suite: 'flis')

     soper_exec = executable('soper_test',
                             'soper_test.c',
                             dependencies : [siril_dep, criterion_dep],
                             link_args : siril_link_arg,
                             c_args : siril_c_flag,
                             cpp_args : siril_cpp_flag)

     test('soper_test', soper_exec, suite: 'arithmetic')

     nde_history_exec = executable('test_nde_history',
                                   'test_nde_history.c',
                                   dependencies : [siril_dep, criterion_dep],
                                   link_args : siril_link_arg,
                                   c_args : siril_c_flag,
                                   cpp_args : siril_cpp_flag)

     test('test_nde_history', nde_history_exec)

     op_descriptor_exec = executable('test_op_descriptor',
                                     'test_op_descriptor.c',
                                     dependencies : [siril_dep, criterion_dep],
                                     link_args : siril_link_arg,
                                     c_args : siril_c_flag,
                                     cpp_args : siril_cpp_flag)

     test('test_op_descriptor', op_descriptor_exec)

     nde_serializers_exec = executable('test_nde_serializers',
                                       'test_nde_serializers.c',
                                       dependencies : [siril_dep, criterion_dep],
                                       link_args : [siril_link_arg, '-Wl,--unresolved-symbols=ignore-all'],
                                       c_args : siril_c_flag,
                                       cpp_args : siril_cpp_flag)

     test('test_nde_serializers', nde_serializers_exec)

     wavelet_denoise_exec = executable('wavelet_denoise_test',
                                       'wavelet_denoise_test.c',
                                       dependencies : [siril_dep, criterion_dep],
                                       link_args : [siril_link_arg, '-Wl,--unresolved-symbols=ignore-all'],
                                       c_args : siril_c_flag,
                                       cpp_args : siril_cpp_flag)

     test('wavelet_denoise_test', wavelet_denoise_exec, suite: 'science')

     photometry_exec = executable('photometry_test',
                                  'photometry_test.c',
                                  dependencies : [siril_dep, criterion_dep],
                                  link_args : siril_link_arg,
                                  c_args : siril_c_flag,
                                  cpp_args : siril_cpp_flag)

     test('photometry_test', photometry_exec, suite: 'science')

     rejection_exec = executable('rejection_test',
                                 'rejection_test.c',
                                 dependencies : [siril_dep, criterion_dep],
                                 link_args : siril_link_arg,
                                 c_args : siril_c_flag,
                                 cpp_args : siril_cpp_flag)

     test('rejection_test', rejection_exec, suite: 'arithmetic')

     ser_exec = executable('ser_test',
                           'ser_test.c',
                           dependencies : [siril_dep, criterion_dep],
                           link_args : siril_link_arg,
                           c_args : siril_c_flag,
                           cpp_args : siril_cpp_flag)

     test('ser_test', ser_exec)

     healpix_cat_exec = executable('healpix_cat_test',
                                   'healpix_cat_test.c',
                                   dependencies : [siril_dep, criterion_dep],
                                   link_args : [siril_link_arg, '-Wl,--unresolved-symbols=ignore-all'],
                                   c_args : siril_c_flag,
                                   cpp_args : siril_cpp_flag)

     test('healpix_cat_test', healpix_cat_exec)

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
                  'mpp_pipeline_test']
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

     homography_exec = executable('homography_test',
                                  'homography_test.c',
                                  dependencies : [siril_dep, criterion_dep],
                                  link_args : [siril_link_arg, '-Wl,--unresolved-symbols=ignore-all'],
                                  c_args : siril_c_flag,
                                  cpp_args : siril_cpp_flag)

     test('homography_test', homography_exec, suite: 'arithmetic')

     gps_exec = executable('qhy_gps_test',
                           ['qhy_gps_test.c', 'download_files.c'],
                           dependencies : [siril_dep, criterion_dep],
                           link_args : siril_link_arg,
                           c_args : siril_c_flag,
                           cpp_args : siril_cpp_flag)

     test('qhy_gps_test', gps_exec)

     streaks_exec = executable('streaks_test',
                           ['streaks_test.c', 'download_files.c'],
                           dependencies : [siril_dep, criterion_dep],
                           link_args : siril_link_arg,
                           c_args : siril_c_flag,
                           cpp_args : siril_cpp_flag)

     test('streaks_test', streaks_exec)

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

     fits_scaling_exec = executable('fits_scaling_test',
                                    'fits_scaling_test.c',
                                    dependencies : [siril_dep, criterion_dep],
                                    link_args : siril_link_arg,
                                    c_args : siril_c_flag,
                                    cpp_args : siril_cpp_flag)

     test('fits_scaling_test', fits_scaling_exec)

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

     test_flis_geometry_ops_exec = executable('test_flis_geometry_ops',
                                'test_flis_geometry_ops.c',
                                dependencies : [siril_dep, criterion_dep],
                                link_args : [siril_link_arg, '-Wl,--unresolved-symbols=ignore-all'],
                                c_args : siril_c_flag,
                                cpp_args : siril_cpp_flag)
     test('test_flis_geometry_ops', test_flis_geometry_ops_exec, suite: 'flis')

     test_flis_mask_route_exec = executable('test_flis_mask_route',
                                'test_flis_mask_route.c',
                                dependencies : [siril_dep, criterion_dep],
                                link_args : [siril_link_arg, '-Wl,--unresolved-symbols=ignore-all'],
                                c_args : siril_c_flag,
                                cpp_args : siril_cpp_flag)
     test('test_flis_mask_route', test_flis_mask_route_exec, suite: 'flis')

     test_flis_layers_match_exec = executable('test_flis_layers_match',
                                'test_flis_layers_match.c',
                                dependencies : [siril_dep, criterion_dep],
                                link_args : [siril_link_arg, '-Wl,--unresolved-symbols=ignore-all'],
                                c_args : siril_c_flag,
                                cpp_args : siril_cpp_flag)
     test('test_flis_layers_match', test_flis_layers_match_exec, suite: 'flis')

     test_flis_register_layers_exec = executable('test_flis_register_layers',
                                'test_flis_register_layers.c',
                                dependencies : [siril_dep, criterion_dep],
                                link_args : [siril_link_arg, '-Wl,--unresolved-symbols=ignore-all'],
                                c_args : siril_c_flag,
                                cpp_args : siril_cpp_flag)
     test('test_flis_register_layers', test_flis_register_layers_exec, suite: 'flis')

     test_flis_register_nde_exec = executable('test_flis_register_nde',
                                'test_flis_register_nde.c',
                                dependencies : [siril_dep, criterion_dep],
                                link_args : [siril_link_arg, '-Wl,--unresolved-symbols=ignore-all'],
                                c_args : siril_c_flag,
                                cpp_args : siril_cpp_flag)
     test('test_flis_register_nde', test_flis_register_nde_exec, suite: 'flis')

     test_flis_canvas_exec = executable('test_flis_canvas',
                                'test_flis_canvas.c',
                                dependencies : [siril_dep, criterion_dep],
                                link_args : [siril_link_arg, '-Wl,--unresolved-symbols=ignore-all'],
                                c_args : siril_c_flag,
                                cpp_args : siril_cpp_flag)
     test('test_flis_canvas', test_flis_canvas_exec, suite: 'flis')

     test_flis_canvas_cmd_exec = executable('test_flis_canvas_cmd',
                                'test_flis_canvas_cmd.c',
                                dependencies : [siril_dep, criterion_dep],
                                link_args : [siril_link_arg, '-Wl,--unresolved-symbols=ignore-all'],
                                c_args : siril_c_flag,
                                cpp_args : siril_cpp_flag)
     test('test_flis_canvas_cmd', test_flis_canvas_cmd_exec, suite: 'flis')

     test_flis_sparse_exec = executable('test_flis_sparse',
                                'test_flis_sparse.c',
                                dependencies : [siril_dep, criterion_dep],
                                link_args : [siril_link_arg, '-Wl,--unresolved-symbols=ignore-all'],
                                c_args : siril_c_flag,
                                cpp_args : siril_cpp_flag)
     test('test_flis_sparse', test_flis_sparse_exec, suite: 'flis')

     test_flis_stress_exec = executable('test_flis_stress',
                                'test_flis_stress.c',
                                dependencies : [siril_dep, criterion_dep],
                                link_args : [siril_link_arg, '-Wl,--unresolved-symbols=ignore-all'],
                                c_args : siril_c_flag,
                                cpp_args : siril_cpp_flag)
     test('test_flis_stress', test_flis_stress_exec, suite: 'flis')

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

     test_nde_snapstore_exec = executable('test_nde_snapstore',
                                'test_nde_snapstore.c',
                                dependencies : [siril_dep, criterion_dep],
                                link_args : [siril_link_arg, '-Wl,--unresolved-symbols=ignore-all'],
                                c_args : siril_c_flag,
                                cpp_args : siril_cpp_flag)
     test('test_nde_snapstore', test_nde_snapstore_exec, suite: 'flis')

     test_nde_graph_exec = executable('test_nde_graph',
                                'test_nde_graph.c',
                                dependencies : [siril_dep, criterion_dep],
                                link_args : [siril_link_arg, '-Wl,--unresolved-symbols=ignore-all'],
                                c_args : siril_c_flag,
                                cpp_args : siril_cpp_flag)
     test('test_nde_graph', test_nde_graph_exec, suite: 'flis')

     test_nde_composite_exec = executable('test_nde_composite',
                                'test_nde_composite.c',
                                dependencies : [siril_dep, criterion_dep],
                                link_args : [siril_link_arg, '-Wl,--unresolved-symbols=ignore-all'],
                                c_args : siril_c_flag,
                                cpp_args : siril_cpp_flag)
     test('test_nde_composite', test_nde_composite_exec, suite: 'flis')

     test_nde_retention_exec = executable('test_nde_retention',
                                'test_nde_retention.c',
                                dependencies : [siril_dep, criterion_dep],
                                link_args : [siril_link_arg, '-Wl,--unresolved-symbols=ignore-all'],
                                c_args : siril_c_flag,
                                cpp_args : siril_cpp_flag)
     test('test_nde_retention', test_nde_retention_exec, suite: 'flis')

     test_nde_compositing_exec = executable('test_nde_compositing',
                                'test_nde_compositing.c',
                                dependencies : [siril_dep, criterion_dep],
                                link_args : [siril_link_arg, '-Wl,--unresolved-symbols=ignore-all'],
                                c_args : siril_c_flag,
                                cpp_args : siril_cpp_flag)
     test('test_nde_compositing', test_nde_compositing_exec, suite: 'flis')

     test_nde_replay_exec = executable('test_nde_replay',
                                'test_nde_replay.c',
                                dependencies : [siril_dep, criterion_dep],
                                link_args : [siril_link_arg, '-Wl,--unresolved-symbols=ignore-all'],
                                c_args : siril_c_flag,
                                cpp_args : siril_cpp_flag)
     test('test_nde_replay', test_nde_replay_exec, suite: 'flis')

     test_nde_conductor_exec = executable('test_nde_conductor',
                                'test_nde_conductor.c',
                                dependencies : [siril_dep, criterion_dep],
                                link_args : [siril_link_arg, '-Wl,--unresolved-symbols=ignore-all'],
                                c_args : siril_c_flag,
                                cpp_args : siril_cpp_flag)
     test('test_nde_conductor', test_nde_conductor_exec, suite: 'flis')

     test_nde_script_scope_exec = executable('test_nde_script_scope',
                                'test_nde_script_scope.c',
                                dependencies : [siril_dep, criterion_dep],
                                link_args : [siril_link_arg, '-Wl,--unresolved-symbols=ignore-all'],
                                c_args : siril_c_flag,
                                cpp_args : siril_cpp_flag)
     test('test_nde_script_scope', test_nde_script_scope_exec, suite: 'flis')

     test_nde_tier_c_exec = executable('test_nde_tier_c',
                                'test_nde_tier_c.c',
                                dependencies : [siril_dep, criterion_dep],
                                link_args : [siril_link_arg, '-Wl,--unresolved-symbols=ignore-all'],
                                c_args : siril_c_flag,
                                cpp_args : siril_cpp_flag)
     test('test_nde_tier_c', test_nde_tier_c_exec, suite: 'flis')

     test_nde_checkpoint_exec = executable('test_nde_checkpoint',
                                'test_nde_checkpoint.c',
                                dependencies : [siril_dep, criterion_dep],
                                link_args : [siril_link_arg, '-Wl,--unresolved-symbols=ignore-all'],
                                c_args : siril_c_flag,
                                cpp_args : siril_cpp_flag)
     test('test_nde_checkpoint', test_nde_checkpoint_exec, suite: 'flis')

     test_nde_persist_exec = executable('test_nde_persist',
                                'test_nde_persist.c',
                                dependencies : [siril_dep, criterion_dep],
                                link_args : [siril_link_arg, '-Wl,--unresolved-symbols=ignore-all'],
                                c_args : siril_c_flag,
                                cpp_args : siril_cpp_flag)
     test('test_nde_persist', test_nde_persist_exec, suite: 'flis')

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

     test_flis_open_exec = executable('test_flis_open',
                                'test_flis_open.c',
                                dependencies : [siril_dep, criterion_dep],
                                link_args : [siril_link_arg, '-Wl,--unresolved-symbols=ignore-all'],
                                c_args : siril_c_flag,
                                cpp_args : siril_cpp_flag)
     test('test_flis_open', test_flis_open_exec, suite: 'flis')

     test_flis_gpu_compose_bake_exec = executable('test_flis_gpu_compose_bake',
                                'test_flis_gpu_compose_bake.c',
                                dependencies : [siril_dep, criterion_dep],
                                link_args : [siril_link_arg, '-Wl,--unresolved-symbols=ignore-all'],
                                c_args : siril_c_flag,
                                cpp_args : siril_cpp_flag)
     test('test_flis_gpu_compose_bake', test_flis_gpu_compose_bake_exec, suite: 'flis')

     test_icc_state_exec = executable('test_icc_state',
                                'test_icc_state.c',
                                dependencies : [siril_dep, criterion_dep],
                                link_args : [siril_link_arg, '-Wl,--unresolved-symbols=ignore-all'],
                                c_args : siril_c_flag,
                                cpp_args : siril_cpp_flag)
     test('test_icc_state', test_icc_state_exec, suite: 'flis')

     # Performance testing
     sorting_perf_exec = executable('sorting_perf',
                                   'measure_sorting.c',
                                   dependencies : siril_dep,
                                   link_args : siril_link_arg,
                                   c_args : siril_c_flag,
                                   cpp_args : siril_cpp_flag)

     test('sorting_perf', sorting_perf_exec, suite: 'perfs', protocol: 'exitcode', is_parallel : false)

endif

