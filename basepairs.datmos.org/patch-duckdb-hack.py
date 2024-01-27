import os, sys, re, pathlib

def patch_duckdb_hack(directory):
    files = list(pathlib.Path(directory).rglob('*.js'))
    done = None
    for file in files:
        content = file.read_text('utf-8')
        if '.headers.get("content-length")' in content and 'new Response' in content and 'WebAssembly.instantiateStreaming' in content:
            if done:
                raise Exception(f'Multiple files found: {done} and {file}')
            done = file
            content = re.sub(r'return new Response\(', 'return duckdb_init_response_content_type_hack(', content)
            if content.count('duckdb_init_response_content_type_hack') != 1:
                raise Exception(f'Failed to patch {file}')
            content = """
                function duckdb_init_response_content_type_hack(...args) {
                    const response = new Response(...args)
                    response.headers.set("content-type", "application/wasm")
                    return response
                }\n""" + content
            print(f'Patched DuckDB {file}')
            
            file.write_text(content, 'utf-8')

    if not done:
        print('No DuckDB files found to patch')

if __name__ == '__main__':
    patch_duckdb_hack(sys.argv[1])

