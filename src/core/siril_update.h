/*
 * This file is part of Siril, an astronomy image processor.
 * Copyright (C) 2005-2011 Francois Meyer (dulle at free.fr)
 * Copyright (C) 2012-2025 team free-astro (see more in AUTHORS file)
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
#ifndef SRC_CORE_SIRIL_UPDATE_H_
#define SRC_CORE_SIRIL_UPDATE_H_

int compare_version(version_number v1, version_number v2);
version_number get_version_number_from_string(const gchar *string);
version_number get_current_version_number();

#if ( defined(HAVE_LIBCURL)  || defined(HAVE_LIBGIT2) )

void siril_check_updates(gboolean verbose);
void siril_check_spcc_mirrors(gboolean verbose);
void siril_check_notifications(gboolean verbose);

#endif

#endif /* SRC_CORE_SIRIL_UPDATE_H_ */
