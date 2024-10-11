import re

# Function to extract API info and output as RST
def extract_api_to_rst(c_file_path, rst_output_path):
    # Regex patterns to match methods, getters/setters, and type definitions
    method_pattern = re.compile(r'{"(\w+)",\s+\(PyCFunction\)(\w+),\s+(METH_[\w|]+),\s+N_\("(.+?)"\)\s*},?')
    type_pattern = re.compile(r'PyTypeObject\s+(\w+)\s+=\s+{.*?\.tp_name\s+=\s+"(.+?)".*?\.tp_doc\s+=\s+N_\("(.+?)"\).*?\.tp_methods\s+=\s+(\w+).*?\.tp_getset\s+=\s+(\w+).*?}', re.DOTALL)

    # Extracted data
    api_documentation = []

    with open(c_file_path, 'r') as c_file:
        content = c_file.read()

        # Extract methods from PyMethodDef
        methods = method_pattern.findall(content)
        if methods:
            api_documentation.append("**Methods**\n")
            for method in methods:
                name, func, flags, doc = method
                api_documentation.append(f".. function:: {name}()\n")
                api_documentation.append(f"   :flags: {flags}\n")
                api_documentation.append(f"   {doc}\n\n")

        # Extract types from PyTypeObject
        types = type_pattern.findall(content)
        if types:
            api_documentation.append("**Types**\n")
            for type_obj in types:
                type_name, tp_name, tp_doc, tp_methods, tp_getset = type_obj
                api_documentation.append(f".. class:: {tp_name}\n")
                api_documentation.append(f"   {tp_doc}\n\n")

                # Look for methods associated with the type
                method_block_pattern = re.compile(rf"static PyMethodDef {tp_methods}\[\].*?\{{(.*?)\}};", re.DOTALL)
                method_block = method_block_pattern.search(content)
                if method_block:
                    api_documentation.append(f"   **Methods for {tp_name}:**\n")
                    methods = method_pattern.findall(method_block.group(1))
                    for method in methods:
                        name, func, flags, doc = method
                        api_documentation.append(f"   - **{name}**: {doc}\n")

                # Look for getters/setters associated with the type
                getset_block_pattern = re.compile(rf"static PyGetSetDef {tp_getset}\[\].*?\{{(.*?)\}};", re.DOTALL)
                getset_block = getset_block_pattern.search(content)
                if getset_block:
                    api_documentation.append(f"   **Getters and Setters for {tp_name}:**\n")
                    getset_pattern = re.compile(r'{"(\w+)",\s+\(getter\)(\w+),\s+\(setter\)(\w+),\s+N_\("(.+?)"\)\s*},?')
                    getters_setters = getset_pattern.findall(getset_block.group(1))
                    for getter_setter in getters_setters:
                        name, getter, setter, doc = getter_setter
                        api_documentation.append(f"   - **{name}**: {doc} (getter: {getter}, setter: {setter})\n")

    # Write the extracted data to an RST file
    with open(rst_output_path, 'w') as rst_file:
        rst_file.write("\n".join(api_documentation))

    print(f"API documentation written to {rst_output_path}")

# Example usage
c_file = "siril_python.c"
rst_output = "siril_python_api.rst"
extract_api_to_rst(c_file, rst_output)
