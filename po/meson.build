################################################################################
# Meson translation file
################################################################################

################################################################################
# Use i18n module

po_dir = meson.current_source_dir()
i18n.gettext(meson.project_name(),
  preset: 'glib',
  # Silence xgettext's "embedded URL" warnings: the remaining URLs are
  # intentional <a href> markup links (help text, GtkLabel markup) that
  # cannot be moved out of the translatable string.
  args: ['--no-check=url'],
)