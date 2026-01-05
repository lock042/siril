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
#ifndef SRC_GUI_REGISTRATION_H_
#define SRC_GUI_REGISTRATION_H_

typedef enum {
	REG_PAGE_GLOBAL,
	REG_PAGE_COMET,
	REG_PAGE_3_STARS,
	REG_PAGE_KOMBAT,
	REG_PAGE_MISC
} reg_notebook_page;

// 3 stars GUI
void reset_3stars();
int _3stars_get_number_selected_stars();
gboolean _3stars_check_selection();

// General GUI
int get_registration_layer_from_GUI(const sequence *seq);
void update_reg_interface(gboolean dont_change_reg_radio);
gboolean end_register_idle(gpointer p);
void initialize_registration_methods();

#endif /* SRC_GUI_REGISTRATION_H_ */
