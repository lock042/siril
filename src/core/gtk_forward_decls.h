/*
 * This file is part of Siril, an astronomy image processor.
 * Copyright (C) 2005-2011 Francois Meyer (dulle at free.fr)
 * Copyright (C) 2012-2026 team free-astro (see more in AUTHORS file)
 * Reference site is https://siril.org
 *
 * Siril is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * Siril is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Siril. If not, see <http://www.gnu.org/licenses/>.
 */

/*
 * Minimal forward declarations for GTK/GDK types that appear as
 * pointer-only parameters or fields in non-GUI headers.
 *
 * Non-GUI code may name these types and pass opaque pointers through,
 * but must never dereference them.  These declarations are skipped when
 * the full <gtk/gtk.h> has already been included (the __GTK_H__ guard).
 *
 * Introduced in plan step 1.4 to allow non-GUI headers to compile after
 * <gtk/gtk.h> was removed from core/siril.h.  Each use site is a
 * temporary violation that will be eliminated in the phases listed below.
 */

#ifndef SIRIL_GTK_FORWARD_DECLS_H
#define SIRIL_GTK_FORWARD_DECLS_H

#ifndef __GTK_H__

/* GdkPixbuf — used as return/parameter type in several headers.
 * ser.c thumbnail use moved to gui/dialog_preview.c; others (proto.h,
 * siril_plot.h, image_format_fits.h) still need this forward declaration. */
typedef struct _GdkPixbuf GdkPixbuf;

/* GtkEntry — parameter of search_object() in gui/search_objects.h. */
typedef struct _GtkEntry GtkEntry;

/* GtkButton — GTK callback declared in algos/astrometry_solver.h.
 * TODO (plan step 3): move callback declaration to gui/astrometry_solver.c. */
typedef struct _GtkButton GtkButton;

/* GtkFileChooser — wrapper helpers declared in core/proto.h.
 * TODO (plan 4+): move file-chooser helpers into gui/ wrappers.             */
typedef struct _GtkFileChooser GtkFileChooser;

#endif /* __GTK_H__ */

#endif /* SIRIL_GTK_FORWARD_DECLS_H */
