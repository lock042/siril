#ifndef _OP_DESCRIPTORS_H_
#define _OP_DESCRIPTORS_H_

/**
 * \file op_descriptors.h
 * \brief Extern declarations for every op_descriptor in the codebase.
 *
 * The declarations are generated from the single canonical list in
 * op_descriptors.def via the X-macro idiom, so this header and the registry
 * array in op_descriptor.c cannot drift out of sync.  Include this header
 * wherever a descriptor is referenced (a construction site \c args->op =
 * &op_desc_x, or the module that defines it).
 *
 * The descriptor DEFINITIONS still live next to their hooks in the owning
 * modules; only the declarations are centralised here.
 */

#include "core/op_descriptor.h"   /* the op_descriptor type */

#ifdef __cplusplus
extern "C" {
#endif

#define OP_DESC(name) extern const op_descriptor name;
#include "core/op_descriptors.def"
#undef OP_DESC

#ifdef __cplusplus
}
#endif

#endif /* _OP_DESCRIPTORS_H_ */
