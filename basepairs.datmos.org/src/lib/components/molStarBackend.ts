import { renderReact18 } from 'molstar/lib/mol-plugin-ui/react18'
import { createPluginUI } from 'molstar/lib/mol-plugin-ui'


export async function createComponent(element: HTMLElement) {

    const ui = createPluginUI({
        target: element,
        render: renderReact18,
        spec: {
            
        }
    })
}
