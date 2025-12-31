/**
 * Utility functions for exporting data and visualizations
 */

/**
 * Export data to CSV format
 */
export const exportToCSV = (
  data: Array<Record<string, any>>,
  filename: string = 'export.csv'
): void => {
  if (data.length === 0) {
    alert('No data to export');
    return;
  }

  const headers = Object.keys(data[0]);
  const csvContent = [
    headers.join(','),
    ...data.map(row =>
      headers.map(header => {
        const value = row[header];
        // Escape quotes and wrap in quotes if contains comma, newline, or quote
        if (value === null || value === undefined) return '""';
        const stringValue = String(value).replace(/"/g, '""');
        if (stringValue.includes(',') || stringValue.includes('\n') || stringValue.includes('"')) {
          return `"${stringValue}"`;
        }
        return stringValue;
      }).join(',')
    )
  ].join('\n');

  const blob = new Blob([csvContent], { type: 'text/csv;charset=utf-8;' });
  const link = document.createElement('a');
  const url = URL.createObjectURL(blob);
  link.setAttribute('href', url);
  link.setAttribute('download', filename);
  link.style.visibility = 'hidden';
  document.body.appendChild(link);
  link.click();
  document.body.removeChild(link);
  URL.revokeObjectURL(url);
};

/**
 * Export data to JSON format
 */
export const exportToJSON = (
  data: any,
  filename: string = 'export.json'
): void => {
  const jsonContent = JSON.stringify(data, null, 2);
  const blob = new Blob([jsonContent], { type: 'application/json;charset=utf-8;' });
  const link = document.createElement('a');
  const url = URL.createObjectURL(blob);
  link.setAttribute('href', url);
  link.setAttribute('download', filename);
  link.style.visibility = 'hidden';
  document.body.appendChild(link);
  link.click();
  document.body.removeChild(link);
  URL.revokeObjectURL(url);
};

/**
 * Export sequences to FASTA format
 */
export const exportToFASTA = (
  sequences: Array<{ name: string; sequence: string }>,
  filename: string = 'sequences.fasta'
): void => {
  if (sequences.length === 0) {
    alert('No sequences to export');
    return;
  }

  const fastaContent = sequences
    .map(seq => `>${seq.name}\n${seq.sequence}`)
    .join('\n');

  const blob = new Blob([fastaContent], { type: 'text/plain;charset=utf-8;' });
  const link = document.createElement('a');
  const url = URL.createObjectURL(blob);
  link.setAttribute('href', url);
  link.setAttribute('download', filename);
  link.style.visibility = 'hidden';
  document.body.appendChild(link);
  link.click();
  document.body.removeChild(link);
  URL.revokeObjectURL(url);
};

/**
 * Type definition for html2canvas function
 */
type Html2Canvas = (element: HTMLElement, options?: any) => Promise<HTMLCanvasElement>;

/**
 * Export element as PNG image
 */
export const exportElementAsPNG = async (
  element: HTMLElement,
  filename: string = 'export.png',
  options?: { backgroundColor?: string; scale?: number }
): Promise<void> => {
  try {
    // Dynamic import to avoid loading html2canvas if not needed
    // Using type assertion since html2canvas types are available but TypeScript may not resolve them correctly
    const html2canvasModule = await import('html2canvas' as any);
    const html2canvas: Html2Canvas = html2canvasModule.default || html2canvasModule;
    
    const canvas = await html2canvas(element, {
      backgroundColor: options?.backgroundColor || '#ffffff',
      scale: options?.scale || 2,
      logging: false,
      useCORS: true,
    });

    canvas.toBlob((blob: Blob | null) => {
      if (!blob) {
        alert('Failed to export image');
        return;
      }
      const url = URL.createObjectURL(blob);
      const link = document.createElement('a');
      link.setAttribute('href', url);
      link.setAttribute('download', filename);
      link.style.visibility = 'hidden';
      document.body.appendChild(link);
      link.click();
      document.body.removeChild(link);
      URL.revokeObjectURL(url);
    }, 'image/png');
  } catch (error) {
    console.error('Error exporting PNG:', error);
    alert('Failed to export image. Please ensure html2canvas is installed. Run: npm install html2canvas');
  }
};

/**
 * Export element as SVG
 */
export const exportElementAsSVG = (
  element: HTMLElement,
  filename: string = 'export.svg'
): void => {
  // Clone the element to avoid modifying the original
  const clonedElement = element.cloneNode(true) as HTMLElement;
  
  // Get computed styles and apply inline styles
  const computedStyles = window.getComputedStyle(element);
  const style = document.createElement('style');
  style.textContent = `
    * {
      font-family: ${computedStyles.fontFamily};
      font-size: ${computedStyles.fontSize};
    }
  `;
  
  const svg = document.createElementNS('http://www.w3.org/2000/svg', 'svg');
  svg.setAttribute('xmlns', 'http://www.w3.org/2000/svg');
  svg.setAttribute('width', element.offsetWidth.toString());
  svg.setAttribute('height', element.offsetHeight.toString());
  
  const foreignObject = document.createElementNS('http://www.w3.org/2000/svg', 'foreignObject');
  foreignObject.setAttribute('width', '100%');
  foreignObject.setAttribute('height', '100%');
  
  foreignObject.appendChild(style);
  foreignObject.appendChild(clonedElement);
  svg.appendChild(foreignObject);
  
  const serializer = new XMLSerializer();
  const svgString = serializer.serializeToString(svg);
  const blob = new Blob([svgString], { type: 'image/svg+xml;charset=utf-8' });
  const url = URL.createObjectURL(blob);
  const link = document.createElement('a');
  link.setAttribute('href', url);
  link.setAttribute('download', filename);
  link.style.visibility = 'hidden';
  document.body.appendChild(link);
  link.click();
  document.body.removeChild(link);
  URL.revokeObjectURL(url);
};

/**
 * Generate filename with timestamp
 */
export const generateFilename = (
  baseName: string,
  extension: string,
  includeTimestamp: boolean = true
): string => {
  if (includeTimestamp) {
    const timestamp = new Date().toISOString().replace(/[:.]/g, '-').slice(0, -5);
    return `${baseName}_${timestamp}.${extension}`;
  }
  return `${baseName}.${extension}`;
};


