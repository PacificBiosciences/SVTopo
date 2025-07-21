
dc.config.defaultColors(d3.schemeSet1)

// plot constraints
const plotw = 585
const ploth = 150

// table filters
var searchInput;
var variantInput;
var nvariantsChart;
var sizeChart;
var chromChart;
var overlapsChart;
// shows filter impact in modal header
var variantCount;

// used to access filtered table data
var chromDimension
// datatables obj
var variant_table
// crossfilter obj
var ndx
// dimensions
var nvariantsDimension
var searchDimension
var variantDimension

// Add this with other global variables at the top
var isUniqueMode = false;

let currentViewer = null;
function handleKeyPress(event) {
    if (!currentViewer) return;

    // Prevent default scrolling behavior for left/right keys
    if (event.key === 'ArrowLeft' || event.key === 'ArrowRight') {
        event.preventDefault();

        const current = $('tr.selected');
        if (!current.length) return;  // No row is selected

        // Get the current row index
        const currentRowIndex = variant_table.row(current).index();
        if (currentRowIndex === undefined) return;

        try {
            // Get visible rows to ensure we navigate within filtered results
            const visibleRows = variant_table.rows({ search: 'applied' });
            const visibleIndexes = visibleRows.indexes();
            const visibleIndexArray = Array.from(visibleIndexes);
            const currentVisibleIndex = visibleIndexArray.indexOf(currentRowIndex);

            if (event.key === 'ArrowLeft' && currentVisibleIndex > 0) {
                const prevVisibleIndex = visibleIndexArray[currentVisibleIndex - 1];
                if (prevVisibleIndex !== undefined) {
                    // Store reference to current viewer before destroying
                    const oldViewer = currentViewer;
                    currentViewer = null;
                    
                    try {
                        oldViewer.destroy();
                        const success = table_click(prevVisibleIndex, variant_table);
                        if (!success) {
                            // Restore viewer state if navigation failed
                            currentViewer = oldViewer;
                        }
                    } catch (error) {
                        console.error("Error navigating to previous row:", error);
                        // Restore viewer state if navigation failed
                        currentViewer = oldViewer;
                    }
                }
            } else if (event.key === 'ArrowRight' && currentVisibleIndex <= visibleIndexArray.length - 2) {
                const nextVisibleIndex = visibleIndexArray[currentVisibleIndex + 1];
                if (nextVisibleIndex !== undefined) {
                    // Store reference to current viewer before destroying
                    const oldViewer = currentViewer;
                    currentViewer = null;
                    
                    try {
                        oldViewer.destroy();
                        const success = table_click(nextVisibleIndex, variant_table);
                        if (!success) {
                            // Restore viewer state if navigation failed
                            currentViewer = oldViewer;
                        }
                    } catch (error) {
                        console.error("Error navigating to next row:", error);
                        // Restore viewer state if navigation failed
                        currentViewer = oldViewer;
                    }
                }
            }
        } catch (error) {
            console.error("Error handling key navigation:", error);
            // Try to clean up if there was an error
            if (currentViewer) {
                try {
                    currentViewer.destroy();
                } catch (e) { /* ignore error during cleanup */ }
                currentViewer = null;
            }
            // Remove the viewer-open class
            document.body.classList.remove('viewer-open');
        }
    }
}

$('#filter-modal').on('hidden.bs.modal', function () {
    update_table()
})

// Charts are rendered during initialization (like paraViewer)
let chartsRendered = true;

// No-op: charts are rendered during initialization, but keep for backward compatibility
function ensureChartsRendered() {
    return chartsRendered;
}

function showImageFromHash() {
    if (window.location.hash) {
        try {
            // Format is #chrom:start-end:sample
            const hash = window.location.hash.substring(1);

            // More permissive pattern matching to handle potential encoding issues
            const match = hash.match(/^([^:]+):(\d+)-(\d+):(.+)$/);

            if (match) {
                const [_, chrom, startStr, endStr, sample] = match;
                const start = parseInt(startStr);
                const end = parseInt(endStr);


                if (isNaN(start) || isNaN(end)) {
                    console.warn(`Invalid coordinates in hash: ${hash}`);
                    history.pushState(null, '', window.location.pathname);
                    return;
                }

                // Find the matching row in visible/filtered data first, then in all data
                let foundRow = false;
                let matchingRowIndex = -1;

                // First try to find in visible rows
                const visibleRows = variant_table.rows({ search: 'applied' });
                visibleRows.every(function (rowIndex) {
                    const data = this.data();

                    if (data.chrom === chrom &&
                        data.start === start &&
                        data.end === end &&
                        data.samples === sample) {
                        matchingRowIndex = rowIndex;
                        foundRow = true;
                        return false; // break the loop
                    }
                    return true;
                });

                // If not found in visible rows, search in all rows
                if (!foundRow) {
                    variant_table.rows().every(function (rowIndex) {
                        const data = this.data();

                        if (data.chrom === chrom &&
                            data.start === start &&
                            data.end === end &&
                            data.samples === sample) {
                            matchingRowIndex = rowIndex;
                            foundRow = true;
                            return false; // break the loop
                        }
                        return true;
                    });
                }

                // Handle case where a matching row was found - use the index directly
                if (foundRow && matchingRowIndex >= 0) {
                    // Use the row index directly - no need to find the DOM node
                    table_click(matchingRowIndex, variant_table);
                } else {
                    console.warn(`No matching row found for ${chrom}:${start}-${end}:${sample}`);
                    // Clear the hash to prevent repeated errors
                    history.pushState(null, '', window.location.pathname);
                }
            } else {
                console.warn(`Hash does not match expected format: ${hash}`);
                history.pushState(null, '', window.location.pathname);
            }
        } catch (error) {
            console.error("Error processing hash:", error);
            history.pushState(null, '', window.location.pathname);
        }
    }
}

const table_click = (selection, table) => {
    // Ensure we have a valid selection before proceeding
    if (selection == null || !table) {
        console.warn("Invalid selection");
        return false;
    }

    let rowIndex;

    // Handle both DOM nodes and direct row indices
    if (typeof selection === 'number') {
        // It's already a row index
        rowIndex = selection;
    } else {
        // It's a DOM node
        rowIndex = table.row(selection).index();
    }

    if (rowIndex === undefined) {
        console.warn("Could not determine row index");
        return false;
    }

    // Validate row index is within bounds of visible rows
    const checkVisibleRows = table.rows({ search: 'applied' });
    const checkVisibleIndexes = checkVisibleRows.indexes();
    
    // Convert to array if it's not already an array
    const visibleIndexArray = Array.from(checkVisibleIndexes);
    
    if (rowIndex < 0 || !visibleIndexArray.includes(rowIndex)) {
        console.warn(`Row index ${rowIndex} is not in visible rows`);
        return false;
    }

    // Get row data before making any DOM changes
    const row = table.row(rowIndex).data();
    if (!row) {
        console.warn("Selected row has no data");
        return false;
    }

    // Check if we have a valid row with the required properties
    if (!row.chrom || !row.start || !row.end || !row.samples || !row.image_name) {
        console.warn("Selected row is missing required data:",
            !row.chrom ? "chrom" : "",
            !row.start ? "start" : "",
            !row.end ? "end" : "",
            !row.samples ? "samples" : "",
            !row.image_name ? "image_name" : "");
        return false;
    }

    // Try to get the DOM node for visual selection
    const rowNode = table.row(rowIndex).node();

    // Only update visual selection if we have a valid DOM node
    if (rowNode) {
        table.$('tr.selected').removeClass('selected');
        $(rowNode).addClass('selected');
    } else {
        console.log("Row node not available for visual selection");
    }

    // Get references to adjacent rows for navigation buttons using visible rows
    const currentRowIndex = rowIndex; // Use the index we already determined
    const navVisibleRows = table.rows({ search: 'applied' });
    const navVisibleIndexes = navVisibleRows.indexes();
    const navVisibleIndexArray = Array.from(navVisibleIndexes);
    const currentVisibleIndex = navVisibleIndexArray.indexOf(currentRowIndex);
    const hasPrevRow = currentVisibleIndex > 0;
    const hasNextRow = currentVisibleIndex < navVisibleIndexArray.length - 1;


    // Update URL hash with row information including sample
    window.history.pushState(null, '',
        `#${row.chrom}:${row.start}-${row.end}:${row.samples}`);

    let img = new Image()
    img.src = `images/${row.image_name}.${plot_type}`

    img.onerror = function () {
        console.error(`Image not found: ${img.src}`);
        alert(`${img.src} not found`);
    }


    let viewer = new Viewer(img, {
        hidden: function () {
            document.removeEventListener('keydown', handleKeyPress);
            currentViewer = null;
            viewer.destroy();

            // Re-enable default browser keyboard behavior when viewer is closed
            document.body.classList.remove('viewer-open');
        },
        title: function () {
            return `chromosome ${row.chrom} at ${row.start}-${row.end}`
        },
        shown: function () {
            // Disable default browser keyboard behavior when viewer is open
            document.body.classList.add('viewer-open');
        },
        toolbar: {
            zoomIn: 4,
            zoomOut: 4,
            oneToOne: 4,
            reset: 4,
            prev: {
                show: hasPrevRow,
                size: "large",
                click: function () {
                    try {
                        if (!hasPrevRow) return;

                        // Get the previous visible row index
                        const visibleRows = table.rows({ search: 'applied' });
                        const visibleIndexes = visibleRows.indexes();
                        const currentVisibleIndex = visibleIndexes.indexOf(currentRowIndex);
                        const prevVisibleIndex = visibleIndexes[currentVisibleIndex - 1];

                        if (prevVisibleIndex !== undefined) {
                            viewer.destroy();
                            currentViewer = null;
                            table_click(prevVisibleIndex, table);
                        }
                    } catch (error) {
                        console.error("Error navigating to previous image:", error);
                        // Try to clean up
                        try {
                            viewer.destroy();
                        } catch (e) { /* ignore error during cleanup */ }
                        currentViewer = null;
                    }
                }
            },
            play: { show: false },
            next: {
                show: hasNextRow,
                size: "large",
                click: function () {
                    try {
                        if (!hasNextRow) return;

                        // Get the next visible row index
                        const visibleRows = table.rows({ search: 'applied' });
                        const visibleIndexes = visibleRows.indexes();
                        const currentVisibleIndex = visibleIndexes.indexOf(currentRowIndex);
                        const nextVisibleIndex = visibleIndexes[currentVisibleIndex + 1];

                        if (nextVisibleIndex !== undefined) {
                            viewer.destroy();
                            currentViewer = null;
                            table_click(nextVisibleIndex, table);
                        }
                    } catch (error) {
                        console.error("Error navigating to next image:", error);
                        // Try to clean up
                        try {
                            viewer.destroy();
                        } catch (e) { /* ignore error during cleanup */ }
                        currentViewer = null;
                    }
                }
            },
            rotateLeft: { show: false },
            rotateRight: { show: false },
            flipHorizontal: { show: false },
            flipVertical: { show: false },
        },
        transition: false,
        navbar: false,
    })

    // Add the event listener when viewer is shown
    currentViewer = viewer;
    document.addEventListener('keydown', handleKeyPress);

    viewer.show()
    
    return true;
}


function build_table(data) {
    // hide the placeholder and show the datatable
    d3.select('#variant-table-placeholder').property("hidden", true)
    d3.select('#variant-table-div').property("hidden", false)

    let cols = [
        { data: 'chrom', title: 'Chrom' },
        { data: 'start', title: 'Start' },
        { data: 'end', title: 'End' },
        { data: 'svlength', title: 'Size' },
        { data: 'samples', title: 'Sample' },
        { data: 'variant_ids', title: 'Variant IDs' },
        { data: 'nvariants', title: '# of Variants' },
    ]

    variant_table = $("#variant-table").DataTable({
        data: data,
        columns: cols,
        deferRender: true,
        scrollY: '80vh',
        scrollCollapse: true,
        scroller: {
            loadingIndicator: true,
            displayBuffer: 20,
            serverWait: 100
        },
        info: true,
        buttons: [
            {
                extend: 'copyHtml5',
                fieldSeparator: '\t',
                fieldBoundary: '',
                exportOptions: {
                    columns: [0, 1, 2, 3, 4, 5, 6],
                    format: {
                        // Special handling for columns with HTML wrappers
                        body: function (data, row, column, node) {
                            // Get the raw data from the row
                            const rowData = variant_table.row(row).data();

                            // For Sample column (index 4)
                            if (column === 4) {
                                return rowData.samples;
                            }
                            // For variant IDs column (index 5)
                            else if (column === 5) {
                                const variantData = rowData.variant_ids;
                                return Array.isArray(variantData) ? variantData.join(', ') : variantData;
                            }
                            // For # of Variants column (index 6)
                            else if (column === 6) {
                                return rowData.nvariants;
                            }
                            return data;
                        }
                    }
                }
            },
            {
                extend: 'csvHtml5',
                fieldSeparator: '\t',
                fieldBoundary: '',
                extension: '.tsv',
                filename: 'svtopo_variants',
                exportOptions: {
                    columns: [0, 1, 2, 3, 4, 5, 6],
                    format: {
                        // Special handling for columns with HTML wrappers
                        body: function (data, row, column, node) {
                            // Get the raw data from the row
                            const rowData = variant_table.row(row).data();

                            // For Sample column (index 4)
                            if (column === 4) {
                                return rowData.samples;
                            }
                            // For variant IDs column (index 5)
                            else if (column === 5) {
                                const variantData = rowData.variant_ids;
                                return Array.isArray(variantData) ? variantData.join(', ') : variantData;
                            }
                            // For # of Variants column (index 6)
                            else if (column === 6) {
                                return rowData.nvariants;
                            }
                            return data;
                        }
                    }
                }
            }
        ],
        dom: 'Brti',
        // Hide the buttons container but keep it in the DOM
        initComplete: function () {
            $('.dt-buttons').hide();

            // Process hash URL only after table is fully initialized
            setTimeout(function () {
                showImageFromHash();
            }, 200);
        },
        infoCallback: (oSettings, iStart, iEnd, iMax, iTotal, sPre) => {
            return `
        <span class="datatable-info"> 
            <span class="pr-2">Showing <b>${iStart}</b> - <b>${iEnd}</b> of <b>${iTotal}</b> records</span>
            <button type="button" class="btn btn-primary btn-sm" data-toggle="modal" data-target="#filter-modal" title="Show filters">
                <span class="fas fa-filter"></span>
            </button>
            <button type="button" class="btn btn-primary btn-sm mr-2" onclick="resetFilter()" title="Reset variant filter">
                <span class="fas fa-undo"></span>
            </button>
            <span class="dropup">
                <button type="button" class="btn btn-sm btn-primary dropdown-toggle" id="download-menu" title="Save table" data-toggle="dropdown" aria-haspopup="true" aria-expanded="false">
                    <span class="fas fa-save"></span>
                </button>
                <span class="dropdown-menu" aria-labelledby="download-menu">
                    <h6 class="dropdown-header">Save ${iTotal} rows as:</h6>
                    <button class="dropdown-item" type="button" id="tsv-button-download" onclick="tsv_button_click()">
                        TSV
                    </button>
                    <button class="dropdown-item" type="button" id="copy-button-download" onclick="copy_button_click()">
                        Copy
                    </button>
                </span>
            </span>
            <div class="custom-control custom-checkbox d-inline-block ml-2">
                <input type="checkbox" class="custom-control-input" id="uniqueToggle"${isUniqueMode ? ' checked' : ''}>
                <label class="custom-control-label" for="uniqueToggle">Show multi-region entries only once</label>
            </div>
            <div class="custom-control custom-checkbox d-inline-block ml-2">
                <input type="checkbox" class="custom-control-input" id="helpToggle">
                <label class="custom-control-label" for="helpToggle">Show help</label>
            </div>
        </span>
        `
        },
        columnDefs: [
            {
                targets: ([0, 1, 2, 3]),
            },
            {
                targets: 4,  // Sample column
                render: function (data, type, row) {
                    if (type === 'display' && data != null) {
                        return `<a href="#" onclick="filterBySample('${data}', event)" style="color: inherit; text-decoration: none;">${data}</a>`;
                    }
                    return data;
                }
            },
            {
                targets: 5,  // Variant IDs column
                render: function (data, type, row) {
                    if (type === 'display' && data != null) {
                        let variantList = Array.isArray(data) ? data : [data];
                        let maxDisplay = 4;  // Maximum number of variants to display
                        let displayText = variantList.length > maxDisplay ?
                            variantList.slice(0, maxDisplay).join(', ') + '...' :
                            variantList.join(', ');

                        let dropdownItems = variantList.map(id =>
                            `<a class="dropdown-item" onclick="filterByVariantId('${id}', event)">${id}</a>`
                        ).join('');

                        return `<div class="dropdown">
                            <span>${displayText}</span>
                            <div class="dropdown-content">
                                ${dropdownItems}
                            </div>
                        </div>`;
                    }
                    return Array.isArray(data) ? data.join(', ') : data;
                }
            },
            {
                targets: 6,  // # of Variants column
                render: function (data, type, row) {
                    if (type === 'display' && data != null) {
                        return `<a href="#" onclick="filterByNVariants(${data}, event)" style="color: inherit; text-decoration: none;">${data}</a>`;
                    }
                    return data;
                }
            }
        ],
        lengthChange: false,
        order: [[0, 'asc'], [1, 'asc']],
        rowId: function (data) {
            return `row-${data.chrom}-${data.start}-${data.end}-${data.samples.replace(/\s+/g, '-')}`;
        },
    })

    // Move the handler outside build_table, at the same level as other event handlers
    $(document).on('change', '#uniqueToggle', function () {
        isUniqueMode = this.checked;
        // Reset the crossfilter with new data
        ndx.remove();
        ndx.add(isUniqueMode ? unique_data : data);

        // Recreate dimensions using the helper function
        const dims = createDimensions(ndx);
        chromDimension = dims.chrom;
        nvariantsDimension = dims.nvariants;
        searchDimension = dims.search;
        variantDimension = dims.variant;

        // Reset filters and update everything (charts are always rendered)
        dc.filterAll();
        try {
            dc.renderAll();
        } catch (error) {
            // Ignore chart rendering errors - table will still update
        }

        // Update the table with new data
        update_table();
    });

    // Help toggle handler
    $(document).on('change', '#helpToggle', function () {
        const helpDiv = document.getElementById('help-div');
        if (this.checked) {
            helpDiv.style.display = 'block';
        } else {
            helpDiv.style.display = 'none';
        }
    });

    // register table clicks on sample_column
    variant_table.on('click', 'tr', function () {
        table_click(this, variant_table)
    })

    // Handle context menu (right-click)
    setupContextMenu();
}

function tsv_button_click() {
    // Get all data that passes the current filters
    let filteredData = [];
    variant_table.rows({ search: 'applied' }).every(function () {
        const rowData = this.data();
        filteredData.push([
            rowData.chrom,
            rowData.start,
            rowData.end,
            rowData.svlength,
            rowData.samples,
            Array.isArray(rowData.variant_ids) ? rowData.variant_ids.join(', ') : rowData.variant_ids,
            rowData.nvariants
        ]);
    });

    // Create TSV content
    let tsvContent = "Chrom\tStart\tEnd\tSize\tSample\tVariant IDs\t# of Variants\n";
    filteredData.forEach(row => {
        tsvContent += row.join('\t') + '\n';
    });

    // Create and trigger download
    const blob = new Blob([tsvContent], { type: 'text/tab-separated-values' });
    const url = window.URL.createObjectURL(blob);
    const a = document.createElement('a');
    a.style.display = 'none';
    a.href = url;
    a.download = 'svtopo_variants.tsv';
    document.body.appendChild(a);
    a.click();
    window.URL.revokeObjectURL(url);
    document.body.removeChild(a);
}

function copy_button_click() {
    // Get all data that passes the current filters
    let filteredData = [];
    variant_table.rows({ search: 'applied' }).every(function () {
        const rowData = this.data();
        filteredData.push([
            rowData.chrom,
            rowData.start,
            rowData.end,
            rowData.svlength,
            rowData.samples,
            Array.isArray(rowData.variant_ids) ? rowData.variant_ids.join(', ') : rowData.variant_ids,
            rowData.nvariants
        ]);
    });

    // Create TSV content for clipboard
    let tsvContent = "Chrom\tStart\tEnd\tSize\tSample\tVariant IDs\t# of Variants\n";
    filteredData.forEach(row => {
        tsvContent += row.join('\t') + '\n';
    });

    // Copy to clipboard
    const textarea = document.createElement('textarea');
    textarea.value = tsvContent;
    document.body.appendChild(textarea);
    textarea.select();
    document.execCommand('copy');
    document.body.removeChild(textarea);

    // Show a notification that data was copied
    alert(`Copied ${filteredData.length} rows to clipboard`);
}

function update_table() {
    variant_table.clear()
    variant_table.rows.add(chromDimension.top(Infinity))
    variant_table.draw()
}

function remove_empty_bins(source_group) {
    return {
        all: function () {
            return source_group.all().filter(function (d) {
                return d.value != 0
            })
        }
    }
}

// https://jsfiddle.net/gordonwoodhull/g34Ldwaz/8/
// https://github.com/dc-js/dc.js/issues/348
function index_group(group) {
    return {
        all: function () {
            return group.all().map(function (kv, i) {
                return { key: i, value: kv.value }
            })
        }
    }
}

// Function to create all dimensions
function createDimensions(ndx) {
    return {
        chrom: ndx.dimension(d => d.chrom),
        nvariants: ndx.dimension(d => {
            const round = d.nvariants < 10 ? 1 :
                d.nvariants < 100 ? 10 :
                    d.nvariants < 1000 ? 100 :
                        d.nvariants < 10000 ? 1000 : 10000;
            return Math.round(d.nvariants / round) * round;
        }),
        search: ndx.dimension(d => d.samples),
        variant: ndx.dimension(d => Array.isArray(d.variant_ids) ? d.variant_ids.join(' ') : d.variant_ids)
    };
}

$(document).ready(function () {
    // Wrap in try-catch to handle potential errors during initialization
    try {
        ndx = crossfilter(data)
        var all = ndx.groupAll()

        // Initial dimension creation using the helper function
        const dims = createDimensions(ndx);
        chromDimension = dims.chrom;
        nvariantsDimension = dims.nvariants;
        searchDimension = dims.search;
        variantDimension = dims.variant;

        build_table(chromDimension.top(Infinity))
        var chromGroup = chromDimension.group().reduceCount()
        var nonEmptyChromGroup = remove_empty_bins(chromGroup)

        searchInput = dc.textFilterWidget("#sample-search")
        searchInput
            .dimension(searchDimension)
            .on('renderlet', function () {
                d3.selectAll(".dc-text-filter-input")
                    .classed("form-control", true)
                d3.selectAll("#sample-search.dc-chart")
                    .classed("col-12", true)
            })

        variantInput = dc.textFilterWidget("#variant-search")
        variantInput
            .dimension(variantDimension)
            .on('renderlet', function () {
                d3.selectAll(".dc-text-filter-input")
                    .classed("form-control", true)
                d3.selectAll("#variant-search.dc-chart")
                    .classed("col-12", true)
            })

        nvariantsChart = dc.barChart("#nvariants-chart")
        var sizeDimension = ndx.dimension(function (d) {
            var round
            if (d.svlength < 100) {
                round = 100
            } else if (d.svlength < 1000) {
                round = 100
            } else if (d.svlength < 10000) {
                round = 1000
            } else if (d.svlength < 100000) {
                round = 10000
            } else if (d.svlength < 1000000) {
                round = 100000
            } else if (d.svlength < 10000000) {
                round = 1000000
            } else {
                round = 10000000
            }
            return Math.round(d.svlength / round) * round
        })
        var sizeGroup = sizeDimension.group().reduceCount()
        var nonEmptySizeGroup = remove_empty_bins(sizeGroup)
        // for brushing, need to track keys at numeric indexes
        var sizeKeys = nonEmptySizeGroup.all().map(dc.pluck('key')).slice()
        var nvariantsGroup = nvariantsDimension.group().reduceCount()
        var nonEmptyNvariantsGroup = remove_empty_bins(nvariantsGroup)
        var nvariantsKeys = nonEmptyNvariantsGroup.all().map(dc.pluck('key')).slice()

        // number of variants
        nvariantsChart
            .width(plotw).height(ploth).gap(1)
            .margins({ top: 10, right: 50, bottom: 30, left: 40 })
            .x(d3.scaleLinear().domain([0, nvariantsKeys.length]))
            .round(Math.floor)
            .brushOn(true)
            .elasticX(true)
            .dimension(nvariantsDimension)
            .group(index_group(nonEmptyNvariantsGroup))
            .elasticY(true)
            .yAxisLabel('Count')
            .filterPrinter(function (filters) {
                var filter = filters[0]
                return nvariantsKeys[filter[0]] + ' - ' + nvariantsKeys[filter[1]]
            })
        // limit the number of labels along x-axis
        nvariantsChart.xAxis().ticks(20)
        nvariantsChart.yAxis().ticks(5)
        // update labels from keys
        nvariantsChart.xAxis().tickFormat(function (v) {
            return nvariantsKeys[v]
        })
        nvariantsChart.filterHandler(function (dimension, filters) {
            if (filters.length === 0) {
                // the empty case (no filtering)
                dimension.filter(null);
            } else {
                const range = [nvariantsKeys[filters[0][0]], nvariantsKeys[filters[0][1]]];
                dimension.filterRange(range);
            }
            return filters;
        });

        // SV length
        sizeChart = dc.barChart("#size-chart")
        sizeChart
            .width(plotw).height(ploth).gap(1)
            .margins({ top: 10, right: 50, bottom: 30, left: 40 })
            .x(d3.scaleLinear().domain([0, sizeKeys.length]))
            .round(Math.floor)
            .brushOn(true)
            .elasticX(true)
            .dimension(sizeDimension)
            .group(index_group(nonEmptySizeGroup))
            .elasticY(true)
            .yAxisLabel('Count')
            .filterPrinter(function (filters) {
                var filter = filters[0]
                return sizeKeys[filter[0]] + ' - ' + sizeKeys[filter[1]]
            })
            // adds left padding to plots inside filtering panel
            .on('renderlet', function () {
                d3.selectAll("svg")
                    .classed("pl-3", true)
            })
        // limit the number of labels along x-axis
        sizeChart.xAxis().ticks(10)
        sizeChart.yAxis().ticks(5)
        // update labels from keys
        sizeChart.xAxis().tickFormat(function (v) {
            return sizeKeys[v]
        })
        // update the status format for this chart
        sizeChart.filterHandler(function (dimension, filters) {
            if (filters.length === 0) {
                // the empty case (no filtering)
                dimension.filter(null)
            } else {
                dimension.filterRange([sizeKeys[filters[0][0]], sizeKeys[filters[0][1]]])
            }
            return filters
        })

        // chromosome
        chromChart = dc.barChart("#chrom-chart")
        chromChart
            .width(plotw).height(ploth).gap(1)
            .margins({ top: 10, right: 50, bottom: 70, left: 40 })
            .x(d3.scaleBand())
            .xUnits(dc.units.ordinal)
            .yAxisLabel('Count')
            .elasticX(true)
            .elasticY(true)
            .dimension(chromDimension)
            .group(nonEmptyChromGroup)
            .ordering((d) => {
                // Strip out 'chr' prefix if it exists
                let key = d.key.toLowerCase().replace('chr', '');
                let v = parseInt(key);
                if (v) {
                    return v;
                } else {
                    // Handle special cases like 'X', 'Y', 'M'
                    if (key === 'x') return 23;
                    if (key === 'y') return 24;
                    if (key === 'm' || key === 'mt') return 25;
                    return key;
                }
            })
            .on('renderlet', function (chart) {
                chart.selectAll('g.x text')
                    .attr('transform', 'rotate(90)')
                    .attr('text-anchor', 'start')
                    .attr('y', -5)
                    .attr('x', 10);
            });
        chromChart.yAxis().ticks(5)

        variantCount = dc.dataCount("#variant-count");
        variantCount.crossfilter(ndx);
        variantCount.groupAll(ndx.groupAll());
            
        variantCount
            .html({
                some: '<strong>%filter-count</strong> selected out of <strong>%total-count</strong> records' +
                    ' | <a href=\'javascript:resetFilter()\'>Reset All</a>',
                all: '<strong>%total-count</strong> records'
            });

            // Render all components during initialization (like paraViewer)
            dc.renderAll();
            
            // Note: Handled by table's initComplete function
            // No need to call showImageFromHash() here
    } catch (error) {
        console.error("Error initializing application:", error);
    }
})

// Charts are already rendered during initialization

function filterByVariantId(variantId, event) {
    event.stopPropagation();

    // Add custom search function
    $.fn.dataTable.ext.search.push(
        function (settings, data, dataIndex) {
            let variantIds = data[5].split(', '); // Column 5 is Variant IDs (after reordering)
            return variantIds.includes(variantId);
        }
    );

    variant_table.draw();

    // Remove the custom search function
    $.fn.dataTable.ext.search.pop();
}

function resetFilter() {
    // Reset DC.js filters (charts are always rendered)
    dc.filterAll();

    // Clear any DataTables search filters
    variant_table.search('').columns().search('');

    // Clear any custom search functions that might have been added
    while ($.fn.dataTable.ext.search.length > 0) {
        $.fn.dataTable.ext.search.pop();
    }

    // Reset DataTable and update with unfiltered data
    variant_table.clear();

    // Use the appropriate data source based on unique mode
    const sourceData = isUniqueMode && typeof unique_data !== 'undefined' ? unique_data : data;

    // Reset the crossfilter with the appropriate data
    ndx.remove();
    ndx.add(sourceData);

    // Recreate dimensions
    const dims = createDimensions(ndx);
    chromDimension = dims.chrom;
    nvariantsDimension = dims.nvariants;
    searchDimension = dims.search;
    variantDimension = dims.variant;

    // Redraw all DC charts to ensure filters are visually reset
    dc.redrawAll();

    // Update the table with the appropriate data
    variant_table.rows.add(chromDimension.top(Infinity));
    variant_table.draw();

    // Clear the URL hash without triggering a page reload
    history.pushState(null, '', window.location.pathname);

    // Reset the sample search input field
    let searchInput = document.querySelector("#sample-search .dc-text-filter-input");
    if (searchInput) searchInput.value = '';

    // Reset the variant search input field
    let variantInput = document.querySelector("#variant-search .dc-text-filter-input");
    if (variantInput) variantInput.value = '';
}

// Handle browser back/forward buttons
window.addEventListener('popstate', function (event) {
    showImageFromHash();
});

function filterBySample(sample, event) {
    event.preventDefault();
    event.stopPropagation();

    // Add custom search function for exact sample match
    $.fn.dataTable.ext.search.push(
        function (settings, data, dataIndex) {
            return data[4] === sample; // Column 4 is Sample
        }
    );

    // Update the sample search input for visual feedback
    let searchInput = document.querySelector("#sample-search .dc-text-filter-input");
    if (searchInput) searchInput.value = sample;

    // Apply the filter
    variant_table.draw();
}

function filterByNVariants(nvariants, event) {
    event.preventDefault();
    event.stopPropagation();

    // Ensure charts are rendered before trying to filter them
    ensureChartsRendered();

    // Apply the filter for exact match
    nvariantsDimension.filter(nvariants);

    // Redraw all charts to ensure UI consistency
    // Only if charts are rendered
    if (chartsRendered) {
        dc.redrawAll();
    }

    // Update the table with filtered data
    update_table();
}

// Setup context menu functionality
function setupContextMenu() {
    const contextMenu = document.getElementById('context-menu');
    const copyImageUrlItem = document.getElementById('copy-image-url');
    const copyPageUrlItem = document.getElementById('copy-page-url');
    const viewImageItem = document.getElementById('view-image');

    // Prevent default context menu
    document.querySelector('#variant-table').addEventListener('contextmenu', function (e) {
        e.preventDefault();

        // Get the closest tr element (table row)
        const tr = $(e.target).closest('tr')[0];
        if (!tr) return;

        // Get the data for this row
        const rowData = variant_table.row(tr).data();
        if (!rowData) return;

        // Position the context menu
        contextMenu.style.top = `${e.pageY}px`;
        contextMenu.style.left = `${e.pageX}px`;
        contextMenu.style.display = 'block';

        // Store reference to the row for use in menu item clicks
        contextMenu.dataset.rowId = tr.id;
    });

    // Hide context menu when clicking outside
    document.addEventListener('click', function () {
        contextMenu.style.display = 'none';
    });

    // Copy Image URL menu item
    copyImageUrlItem.addEventListener('click', function () {
        const tr = document.getElementById(contextMenu.dataset.rowId);
        if (!tr) return;

        const rowData = variant_table.row(tr).data();
        if (!rowData || !rowData.image_name) return;

        // Create the image URL
        const imageUrl = `${window.location.origin}${window.location.pathname.substring(0, window.location.pathname.lastIndexOf('/') + 1)}${rowData.image_name}.${plot_type}`;

        // Copy to clipboard
        copyToClipboard(imageUrl);
    });

    // Copy Page URL with hash
    copyPageUrlItem.addEventListener('click', function () {
        const tr = document.getElementById(contextMenu.dataset.rowId);
        if (!tr) return;

        const rowData = variant_table.row(tr).data();
        if (!rowData || !rowData.chrom || !rowData.start || !rowData.end || !rowData.samples) return;

        // Create the page URL with hash
        const hash = `#${rowData.chrom}:${rowData.start}-${rowData.end}:${rowData.samples}`;
        const pageUrl = `${window.location.origin}${window.location.pathname}${hash}`;

        // Copy to clipboard
        copyToClipboard(pageUrl);
    });

    // View image menu item
    viewImageItem.addEventListener('click', function () {
        try {
            const tr = document.getElementById(contextMenu.dataset.rowId);
            if (!tr) {
                console.warn("No table row found with ID:", contextMenu.dataset.rowId);
                return;
            }

            console.log("Context menu: viewing image for row", contextMenu.dataset.rowId);

            // Trigger the same action as clicking the row
            table_click(tr, variant_table);
        } catch (error) {
            console.error("Error handling view image action:", error);
        }
    });
}

// Helper function to copy text to clipboard
function copyToClipboard(text) {
    const textarea = document.createElement('textarea');
    textarea.value = text;
    document.body.appendChild(textarea);
    textarea.select();

    try {
        const successful = document.execCommand('copy');
        if (successful) {
            // Show success message
            const successMsg = document.querySelector('.clipboard-success');
            successMsg.style.display = 'block';
            setTimeout(() => {
                successMsg.style.display = 'none';
            }, 2000);
        }
    } catch (err) {
        console.error('Unable to copy to clipboard', err);
        alert('Failed to copy URL to clipboard');
    }

    document.body.removeChild(textarea);
}

// Safe chart reset functions
function resetNVariantsChart() {
    ensureChartsRendered();
    if (chartsRendered && nvariantsChart) {
        nvariantsChart.filterAll();
        dc.redrawAll();
    }
}

function resetSizeChart() {
    ensureChartsRendered();
    if (chartsRendered && sizeChart) {
        sizeChart.filterAll();
        dc.redrawAll();
    }
}

function resetChromChart() {
    ensureChartsRendered();
    if (chartsRendered && chromChart) {
        chromChart.filterAll();
        dc.redrawAll();
    }
}

function resetOverlapsChart() {
    ensureChartsRendered();
    if (chartsRendered && overlapsChart) {
        overlapsChart.filterAll();
        dc.redrawAll();
    }
}

function resetAllCharts() {
    ensureChartsRendered();
    if (chartsRendered) {
        dc.filterAll();
        dc.renderAll();
    }
}
