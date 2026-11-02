#include "KaryogramViewer.h"
#include "generalUtilities.h"
#include "SDDparserVersion.h"

#include <cstdint>
#include <cstdlib>
#include <ctime>
#include <fstream>
#include <iostream>
#include <string>




namespace sddparser
{

namespace
{
    // Encode PNG bytes as base64 so they can be embedded in the HTML as data: URLs.
    std::string base64Encode(const std::vector<unsigned char>& bytes)
    {
        static const char table[] = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789+/";

        std::string encoded;
        encoded.reserve(((bytes.size() + 2) / 3) * 4);

        size_t i = 0;
        for (; i + 2 < bytes.size(); i += 3)                                    // Every 3 input bytes become 4 output characters.
        {
            const unsigned int chunk = (bytes[i] << 16) | (bytes[i + 1] << 8) | bytes[i + 2];
            encoded += table[(chunk >> 18) & 0x3F];
            encoded += table[(chunk >> 12) & 0x3F];
            encoded += table[(chunk >> 6) & 0x3F];
            encoded += table[chunk & 0x3F];
        }

        const size_t remaining = bytes.size() - i;                             // Pad the last 1 or 2 bytes with '='.
        if (remaining == 1)
        {
            const unsigned int chunk = bytes[i] << 16;
            encoded += table[(chunk >> 18) & 0x3F];
            encoded += table[(chunk >> 12) & 0x3F];
            encoded += "==";
        }
        else if (remaining == 2)
        {
            const unsigned int chunk = (bytes[i] << 16) | (bytes[i + 1] << 8);
            encoded += table[(chunk >> 18) & 0x3F];
            encoded += table[(chunk >> 12) & 0x3F];
            encoded += table[(chunk >> 6) & 0x3F];
            encoded += '=';
        }

        return encoded;
    }

    // Escape text placed inside the HTML (file names can contain characters such as & or <).
    std::string htmlEscape(const std::string& text)
    {
        std::string escaped;
        for (char c : text)
        {
            switch (c)
            {
                case '&':  escaped += "&amp;";  break;
                case '<':  escaped += "&lt;";   break;
                case '>':  escaped += "&gt;";   break;
                case '"':  escaped += "&quot;"; break;
                case '\'': escaped += "&#39;";  break;
                default:   escaped += c;
            }
        }
        return escaped;
    }

    // Write text as a JavaScript string literal (with quotes) for use inside the viewer's <script>.
    // Why: file names and titles come from the input file name, and could contain quotes, backslashes or "</script>",
    // which would otherwise break the script. HTML escaping does not apply inside <script>, so JavaScript escaping is used.
    std::string jsString(const std::string& text)
    {
        std::string quoted = "\"";
        for (char c : text)
        {
            switch (c)
            {
                case '\\': quoted += "\\\\"; break;
                case '"':  quoted += "\\\""; break;
                case '<':  quoted += "\\u003c"; break;
                case '\n': quoted += "\\n"; break;
                case '\r': quoted += "\\r"; break;
                default:   quoted += c;
            }
        }
        return quoted + "\"";
    }

#if !defined(_WIN32) && !defined(__APPLE__)
    // Detect Windows Subsystem for Linux (WSL).
    // Why: a program built inside WSL is a Linux program, so it would normally use xdg-open. But WSL usually has no
    // Linux browser installed, so xdg-open fails silently and nothing appears. In WSL we instead hand the file to the
    // Windows browser. WSL sets WSL_DISTRO_NAME, and its kernel version string mentions "microsoft" as a backup check.
    bool isRunningInWSL()
    {
        if (std::getenv("WSL_DISTRO_NAME") != nullptr)
        {
            return true;
        }

        std::ifstream version("/proc/version");
        std::string line;
        std::getline(version, line);
        return toLower(line).find("microsoft") != std::string::npos;
    }

    // Check that WSL can launch Windows programs.
    // Why: "Windows interop" can be switched off in WSL (common on locked-down work laptops). Then wslpath or
    // explorer.exe are missing, and we must fall back to telling the user where the file is.
    bool wslCanLaunchWindowsPrograms()
    {
        return std::system("command -v wslpath > /dev/null 2>&1 && command -v explorer.exe > /dev/null 2>&1") == 0;
    }
#endif

#if defined(__APPLE__)
    // Check whether this session is connected over SSH (only needed on macOS, Linux checks for a display instead).
    // Why: on a remote server, HPC cluster or container reached over SSH, the browser would have to open on the
    // remote machine, which has no screen. SSH sets these variables in every remote session.
    bool isRemoteSSHsession()
    {
        return std::getenv("SSH_CONNECTION") != nullptr || std::getenv("SSH_TTY") != nullptr;
    }
#endif

#if !defined(_WIN32)
    // Quote a path for a POSIX shell (Linux, macOS, WSL).
    // Why: the path goes into a std::system() command line. Input file names can contain spaces, $, quotes or other
    // characters the shell would otherwise interpret, which could break the command or run something unintended.
    // Inside single quotes the shell treats everything literally; an embedded ' is written as '\'' (close, escaped
    // quote, reopen).
    std::string shellQuote(const std::string& path)
    {
        std::string quoted = "'";
        for (char c : path)
        {
            quoted += (c == '\'') ? std::string("'\\''") : std::string(1, c);
        }
        quoted += "'";
        return quoted;
    }
#endif
}


bool writeKaryogramViewer(
    const std::vector<KaryogramImage>& images,
    const std::string& title,
    const std::filesystem::path& htmlPath)
{
    if (images.empty())
    {
        std::cerr << "ERROR: No karyograms to show in the web viewer.\n";
        return false;
    }

    std::ofstream html(htmlPath, std::ios::binary);

    if (!html)
    {
        std::cerr << "ERROR: Failed to create web viewer file " << htmlPath << "\n";
        return false;
    }

    const std::string safeTitle = htmlEscape(title);

    // ------------------------------------------
    // Page header and styling
    // ------------------------------------------
    // The karyogram sits inside a scrollable "viewport". Zooming changes the image's real width and height (not a CSS
    // scale transform), so the browser redraws the SVG at the new size and lines and text stay sharp at every zoom level.
    html << "<!doctype html>\n<html lang=\"en\">\n<head>\n<meta charset=\"utf-8\">\n"
         << "<meta name=\"viewport\" content=\"width=device-width, initial-scale=1\">\n"
         << "<title>" << safeTitle << " karyograms</title>\n"
         << R"HTML(<style>
html, body { height: 100%; margin: 0; }
body { display: flex; flex-direction: column; font-family: sans-serif; background: #e9e9e9; color: #222; }
#titlebar { display: grid; grid-template-columns: 1fr auto 1fr; align-items: center; gap: 12px; padding: 6px 16px 0; background: #fff; }
#version { font-size: 13px; color: #666; white-space: nowrap; }
#titlebar h1 { font-size: 16px; margin: 0; text-align: center; overflow-wrap: anywhere; }
header { display: flex; flex-wrap: wrap; gap: 8px; align-items: center; padding: 6px 16px 8px; background: #fff; border-bottom: 1px solid #ccc; }
.group { display: flex; gap: 4px; align-items: center; }
.spacer { flex: 1; }
button, select { font-size: 14px; padding: 4px 10px; cursor: pointer; }
#counter, #zoomLevel { color: #555; font-size: 14px; min-width: 52px; text-align: center; }
#viewport { flex: 1; overflow: auto; position: relative; cursor: grab; user-select: none; touch-action: none; }
#viewport.dragging { cursor: grabbing; }
#stage { position: relative; }
#frame { position: absolute; background: #fff; box-shadow: 0 1px 6px rgba(0,0,0,0.25); }
#karyogram { display: block; width: 100%; height: 100%; }
.hotspot { position: absolute; cursor: pointer; border-radius: 3px; }
.hotspot:hover { background: rgba(30,110,220,0.15); outline: 1px solid rgba(30,110,220,0.6); }
#popup { position: fixed; z-index: 10; min-width: 200px; max-width: 320px; background: #fff; border: 1px solid #bbb; border-radius: 6px; box-shadow: 0 4px 16px rgba(0,0,0,0.25); padding: 10px 14px; font-size: 14px; }
#popup[hidden] { display: none; }
#popup h2 { font-size: 15px; margin: 0 24px 6px 0; }
#popup h3 { font-size: 13px; margin: 10px 0 4px; color: #444; }
#popup table { border-collapse: collapse; width: 100%; }
#popup td { padding: 2px 0; }
#popup td.count { text-align: right; padding-left: 16px; font-variant-numeric: tabular-nums; }
#popup .none { color: #666; font-style: italic; }
#popup .note { margin-top: 8px; font-size: 11px; color: #777; }
#popupClose { position: absolute; top: 4px; right: 6px; border: none; background: none; font-size: 18px; padding: 0 4px; }
#hint { padding: 4px 16px; font-size: 12px; color: #666; background: #fff; border-top: 1px solid #ccc; }
</style>
</head>
<body>
)HTML";

    // ------------------------------------------
    // Toolbar: switch between karyograms, zoom controls
    // ------------------------------------------
    // Title bar: SDDparser version on the left, input file name in the middle.
    html << "<div id=\"titlebar\"><span id=\"version\">" << htmlEscape(SDDPARSER_VERSION_TEXT) << "</span><h1>" << safeTitle << "</h1><span></span></div>\n"
         << "<header>\n"
         << "<div class=\"group\">\n"
         << "<button id=\"prev\" title=\"Previous karyogram\">&larr; Previous</button>\n"
         << "<select id=\"cell\" title=\"Choose a karyogram\">\n";

    for (size_t i = 0; i < images.size(); i++)
    {
        html << "<option value=\"" << i << "\">" << htmlEscape(images[i].label) << "</option>\n";
    }

    html << R"HTML(</select>
<button id="next" title="Next karyogram">Next &rarr;</button>
<span id="counter"></span>
</div>
<div class="spacer"></div>
<div class="group">
<button id="zoomOut" title="Zoom out (-)">&minus;</button>
<span id="zoomLevel"></span>
<button id="zoomIn" title="Zoom in (+)">+</button>
<button id="fit" title="Fit the whole karyogram in the window (0)">Fit</button>
<button id="actual" title="Show at original size (1)">100%</button>
</div>
<div class="group">
<select id="format" title="File format for downloads">
<option value="png">PNG</option>
<option value="svg">SVG</option>
</select>
<button id="download" title="Download the karyogram shown">Download</button>
<button id="downloadAll" title="Download every karyogram in one .zip file">Download all (.zip)</button>
</div>
</header>
<div id="viewport"><div id="stage"><div id="frame"><img id="karyogram" alt="Karyogram" draggable="false"><div id="hotspots"></div></div></div></div>
<div id="popup" hidden role="dialog"><button id="popupClose" title="Close (Esc)">&times;</button><div id="popupBody"></div></div>
<div id="hint">Click a chromosome label for its summary &middot; scroll or pinch to zoom at the cursor &middot; drag to pan &middot; double-click to zoom in &middot; keys: + &minus; 0 (fit) 1 (100%)</div>
)HTML";

    // ------------------------------------------
    // Embedded karyogram images
    // ------------------------------------------
    // Each image is embedded as a base64 data: URL so the HTML file works on its own, with no other files or server.
    // The SVG is shown and both formats are kept for the download buttons.
    html << "<script>\nconst TITLE = " << jsString(title) << ";\nconst karyograms = [\n";

    for (const KaryogramImage& image : images)
    {
        html << "  { filename: " << jsString(image.filename)
             << ",\n    svg: \"" << (image.svgBytes.empty() ? "" : "data:image/svg+xml;base64," + base64Encode(image.svgBytes))
             << "\",\n    png: \"" << (image.pngBytes.empty() ? "" : "data:image/png;base64," + base64Encode(image.pngBytes))
             << "\",\n    hotspots: [";

        // Clickable chromosome labels: position in image pixels plus the counts shown in the popup.
        for (const KaryogramHotspot& hotspot : image.hotspots)
        {
            html << "\n      { title: " << jsString(hotspot.title)
                 << ", x: " << hotspot.x << ", y: " << hotspot.y
                 << ", w: " << hotspot.width << ", h: " << hotspot.height << ", sections: [";
            for (const KaryogramPopupSection& section : hotspot.sections)
            {
                html << "{ heading: " << jsString(section.heading) << ", counts: [";
                for (const auto& [type, count] : section.counts)
                {
                    html << "[" << jsString(type) << ", " << count << "], ";
                }
                html << "] }, ";
            }
            html << "] },";
        }

        html << " ],\n    popupNote: " << jsString(image.popupNote)
             << ",\n    popupEmptyText: " << jsString(image.popupEmptyText) << " },\n";
    }

    // ------------------------------------------
    // Switching, zoom and pan logic
    // ------------------------------------------
    html << R"HTML(];

const viewport = document.getElementById('viewport');
const stage = document.getElementById('stage');
const frame = document.getElementById('frame');
const img = document.getElementById('karyogram');
const hotspotLayer = document.getElementById('hotspots');
const popup = document.getElementById('popup');
const popupBody = document.getElementById('popupBody');
const select = document.getElementById('cell');
const counter = document.getElementById('counter');
const zoomLevel = document.getElementById('zoomLevel');

const MIN_SCALE = 0.1, MAX_SCALE = 20, STEP = 1.25, MARGIN = 16;
let current = 0;
let scale = 1;
let fitMode = true;            // While true the image refits when the window is resized.

// Scale at which the whole karyogram fits inside the window.
function fitScale() {
  const w = img.naturalWidth || 1, h = img.naturalHeight || 1;
  return Math.min((viewport.clientWidth - 2 * MARGIN) / w, (viewport.clientHeight - 2 * MARGIN) / h);
}

// Offset of the image inside the stage: centred when smaller than the window, otherwise a small margin.
function offsets(s) {
  const w = img.naturalWidth * s, h = img.naturalHeight * s;
  return {
    x: Math.max(MARGIN, (viewport.clientWidth - w) / 2),
    y: Math.max(MARGIN, (viewport.clientHeight - h) / 2),
  };
}

// Apply a new zoom level while keeping the image point under (screenX, screenY) fixed on screen.
function setScale(newScale, screenX, screenY) {
  newScale = Math.min(MAX_SCALE, Math.max(MIN_SCALE, newScale));
  if (screenX === undefined) { screenX = viewport.clientWidth / 2; screenY = viewport.clientHeight / 2; }

  const before = offsets(scale);
  const imageX = (viewport.scrollLeft + screenX - before.x) / scale;    // Point under the cursor, in image pixels.
  const imageY = (viewport.scrollTop + screenY - before.y) / scale;

  scale = newScale;
  const after = offsets(scale);
  const w = img.naturalWidth * scale, h = img.naturalHeight * scale;
  frame.style.width = w + 'px';            // Real layout size, so the SVG is redrawn sharply at this size.
  frame.style.height = h + 'px';           // The clickable labels inside the frame are placed in %, so they follow.
  frame.style.left = after.x + 'px';
  frame.style.top = after.y + 'px';
  stage.style.width = Math.max(viewport.clientWidth, w + 2 * after.x) + 'px';
  stage.style.height = Math.max(viewport.clientHeight, h + 2 * after.y) + 'px';

  viewport.scrollLeft = imageX * scale + after.x - screenX;
  viewport.scrollTop = imageY * scale + after.y - screenY;
  zoomLevel.textContent = Math.round(scale * 100) + '%';
}

function fit() { fitMode = true; setScale(fitScale()); }
function zoomBy(factor, x, y) { fitMode = false; setScale(scale * factor, x, y); }

function show(index) {
  current = (index + karyograms.length) % karyograms.length;
  const first = !img.src;
  closePopup();
  img.onload = () => {
    placeHotspots();
    if (first || fitMode) fit(); else setScale(scale);           // Keep the zoom when switching cells, to compare them.
  };
  img.src = karyograms[current].svg || karyograms[current].png;      // Prefer the sharp vector image.
  select.value = current;
  counter.textContent = (current + 1) + ' of ' + karyograms.length;
}

// Mouse wheel and trackpad pinch (sent as wheel events with ctrlKey) zoom around the cursor.
viewport.addEventListener('wheel', (e) => {
  e.preventDefault();
  const r = viewport.getBoundingClientRect();
  const delta = e.deltaMode === 1 ? e.deltaY * 33 : e.deltaY;           // Line-based scrolling (Firefox) to pixels.
  zoomBy(Math.exp(-delta * (e.ctrlKey ? 0.01 : 0.002)), e.clientX - r.left, e.clientY - r.top);
}, { passive: false });

// Drag with the mouse (or one finger) to pan.
let drag = null;
viewport.addEventListener('pointerdown', (e) => {
  if (e.button !== 0) return;
  drag = { x: e.clientX, y: e.clientY, left: viewport.scrollLeft, top: viewport.scrollTop, moved: false };
  viewport.setPointerCapture(e.pointerId);
  viewport.classList.add('dragging');
});
viewport.addEventListener('pointermove', (e) => {
  if (!drag) return;
  if (Math.abs(e.clientX - drag.x) + Math.abs(e.clientY - drag.y) > 4) drag.moved = true;
  viewport.scrollLeft = drag.left - (e.clientX - drag.x);
  viewport.scrollTop = drag.top - (e.clientY - drag.y);
});
const endDrag = () => { drag = null; viewport.classList.remove('dragging'); };
viewport.addEventListener('pointerup', (e) => {
  // A press and release without dragging is a click. The pointer is captured by the viewport for panning, so the
  // click target is looked up from the position: a chromosome label opens its popup, anywhere else closes it.
  if (drag && !drag.moved) {
    const spot = hotspotAt(e.clientX, e.clientY);
    if (spot) openPopup(spot, e.clientX, e.clientY); else closePopup();
  }
  endDrag();
});
viewport.addEventListener('pointercancel', endDrag);

viewport.addEventListener('dblclick', (e) => {
  if (hotspotAt(e.clientX, e.clientY)) return;                       // Double-clicking a label should not also zoom.
  const r = viewport.getBoundingClientRect();
  zoomBy(2, e.clientX - r.left, e.clientY - r.top);
});

document.getElementById('zoomIn').onclick = () => zoomBy(STEP);
document.getElementById('zoomOut').onclick = () => zoomBy(1 / STEP);
document.getElementById('fit').onclick = fit;
document.getElementById('actual').onclick = () => { fitMode = false; setScale(1); };
document.getElementById('prev').onclick = () => show(current - 1);
document.getElementById('next').onclick = () => show(current + 1);
select.onchange = () => show(Number(select.value));

document.addEventListener('keydown', (e) => {
  if (e.key === 'Escape') { closePopup(); return; }
  if (e.target === select || e.ctrlKey || e.metaKey || e.altKey) return;
  if (e.key === '+' || e.key === '=') zoomBy(STEP);
  else if (e.key === '-' || e.key === '_') zoomBy(1 / STEP);
  else if (e.key === '0') fit();
  else if (e.key === '1') { fitMode = false; setScale(1); }
});

// ------------------------------------------------------------
// Clickable chromosome labels and the summary popup (mutations for SDR, damages for SDD)
// ------------------------------------------------------------
// One invisible box per chromosome label, placed in % of the image size so it stays on the label at every zoom level.
function placeHotspots() {
  hotspotLayer.replaceChildren();
  const w = img.naturalWidth || 1, h = img.naturalHeight || 1;
  karyograms[current].hotspots.forEach((spot, i) => {
    const box = document.createElement('div');
    box.className = 'hotspot';
    box.dataset.index = i;
    box.title = spot.title + ': click for summary';
    box.style.left = (100 * spot.x / w) + '%';
    box.style.top = (100 * spot.y / h) + '%';
    box.style.width = (100 * spot.w / w) + '%';
    box.style.height = (100 * spot.h / h) + '%';
    hotspotLayer.appendChild(box);
  });
}

function hotspotAt(x, y) {
  const box = document.elementFromPoint(x, y)?.closest('.hotspot');
  return box ? karyograms[current].hotspots[Number(box.dataset.index)] : null;
}

function addElement(parent, tag, text, className) {
  const element = document.createElement(tag);
  if (text !== undefined) element.textContent = text;            // textContent, never innerHTML: names are shown as plain text.
  if (className) element.className = className;
  parent.appendChild(element);
  return element;
}

function openPopup(spot, x, y) {
  popupBody.replaceChildren();
  addElement(popupBody, 'h2', spot.title);
  for (const section of spot.sections) {
    if (section.heading) addElement(popupBody, 'h3', section.heading);
    if (section.counts.length === 0) { addElement(popupBody, 'div', karyograms[current].popupEmptyText, 'none'); continue; }
    const table = addElement(popupBody, 'table');
    for (const [type, count] of section.counts) {
      const row = addElement(table, 'tr');
      addElement(row, 'td', type);
      addElement(row, 'td', count, 'count');
    }
  }
  if (karyograms[current].popupNote) addElement(popupBody, 'div', karyograms[current].popupNote, 'note');

  popup.hidden = false;
  // Open next to the click, but keep the whole popup inside the window.
  const r = popup.getBoundingClientRect();
  popup.style.left = Math.max(8, Math.min(x + 12, window.innerWidth - r.width - 8)) + 'px';
  popup.style.top = Math.max(8, Math.min(y + 12, window.innerHeight - r.height - 8)) + 'px';
}

function closePopup() { popup.hidden = true; }
document.getElementById('popupClose').onclick = closePopup;

// ------------------------------------------------------------
// Downloads
// ------------------------------------------------------------
const formatSelect = document.getElementById('format');
for (const option of formatSelect.options) {                       // Only offer formats that were embedded.
  option.disabled = !karyograms.some((k) => k[option.value]);
}
if (formatSelect.selectedOptions[0].disabled) formatSelect.value = 'svg';

// "cell0_karyogram_human.png" -> "cell0_karyogram_human.svg" for SVG downloads.
function downloadName(k, format) {
  return k.filename.replace(/\.(png|svg)$/i, '') + '.' + format;
}

// Decode an embedded base64 data: URL back into the original file bytes.
function dataUrlBytes(dataUrl) {
  const binary = atob(dataUrl.slice(dataUrl.indexOf(',') + 1));
  const bytes = new Uint8Array(binary.length);
  for (let i = 0; i < binary.length; i++) bytes[i] = binary.charCodeAt(i);
  return bytes;
}

// Hand a file to the browser's normal download (it asks where to save, or uses the Downloads folder).
function saveFile(bytes, name, type) {
  const url = URL.createObjectURL(new Blob([bytes], { type }));
  const link = document.createElement('a');
  link.href = url;
  link.download = name;
  document.body.appendChild(link);
  link.click();
  link.remove();
  setTimeout(() => URL.revokeObjectURL(url), 1000);
}

const MIME = { png: 'image/png', svg: 'image/svg+xml' };

document.getElementById('download').onclick = () => {
  const k = karyograms[current], format = formatSelect.value;
  if (!k[format]) { alert('This karyogram has no ' + format.toUpperCase() + ' version.'); return; }
  saveFile(dataUrlBytes(k[format]), downloadName(k, format), MIME[format]);
};

// Build a .zip file in the browser (stored, not compressed: PNGs are already compressed and this keeps the code short).
// Why a zip: browsers block or prompt for many separate downloads at once, a single file avoids that.
const CRC_TABLE = (() => {
  const table = new Uint32Array(256);
  for (let n = 0; n < 256; n++) {
    let c = n;
    for (let k = 0; k < 8; k++) c = c & 1 ? 0xEDB88320 ^ (c >>> 1) : c >>> 1;
    table[n] = c >>> 0;
  }
  return table;
})();

function crc32(bytes) {
  let c = 0xFFFFFFFF;
  for (let i = 0; i < bytes.length; i++) c = CRC_TABLE[(c ^ bytes[i]) & 0xFF] ^ (c >>> 8);
  return (c ^ 0xFFFFFFFF) >>> 0;
}

function makeZip(files) {
  const encoder = new TextEncoder();
  const chunks = [], central = [];
  let offset = 0;
  const now = new Date();
  const dosTime = (now.getHours() << 11) | (now.getMinutes() << 5) | (now.getSeconds() >> 1);
  const dosDate = ((now.getFullYear() - 1980) << 9) | ((now.getMonth() + 1) << 5) | now.getDate();

  for (const file of files) {
    const name = encoder.encode(file.name), crc = crc32(file.bytes), size = file.bytes.length;

    const local = new DataView(new ArrayBuffer(30));
    local.setUint32(0, 0x04034B50, true);   // Local file header signature
    local.setUint16(4, 20, true);           // Version needed to extract
    local.setUint16(6, 0x0800, true);       // Flags: file name is UTF-8
    local.setUint16(8, 0, true);            // Method: stored
    local.setUint16(10, dosTime, true);
    local.setUint16(12, dosDate, true);
    local.setUint32(14, crc, true);
    local.setUint32(18, size, true);        // Compressed size
    local.setUint32(22, size, true);        // Uncompressed size
    local.setUint16(26, name.length, true);
    local.setUint16(28, 0, true);           // Extra field length
    chunks.push(new Uint8Array(local.buffer), name, file.bytes);

    const entry = new DataView(new ArrayBuffer(46));
    entry.setUint32(0, 0x02014B50, true);   // Central directory signature
    entry.setUint16(4, 20, true);           // Version made by
    entry.setUint16(6, 20, true);           // Version needed
    entry.setUint16(8, 0x0800, true);
    entry.setUint16(10, 0, true);
    entry.setUint16(12, dosTime, true);
    entry.setUint16(14, dosDate, true);
    entry.setUint32(16, crc, true);
    entry.setUint32(20, size, true);
    entry.setUint32(24, size, true);
    entry.setUint16(28, name.length, true);
    entry.setUint32(42, offset, true);      // Where this file's local header starts
    central.push(new Uint8Array(entry.buffer), name);

    offset += 30 + name.length + size;
  }

  const centralSize = central.reduce((total, part) => total + part.length, 0);
  const end = new DataView(new ArrayBuffer(22));
  end.setUint32(0, 0x06054B50, true);       // End of central directory signature
  end.setUint16(8, files.length, true);
  end.setUint16(10, files.length, true);
  end.setUint32(12, centralSize, true);
  end.setUint32(16, offset, true);
  return new Blob([...chunks, ...central, new Uint8Array(end.buffer)], { type: 'application/zip' });
}

document.getElementById('downloadAll').onclick = () => {
  const format = formatSelect.value;
  const files = karyograms.filter((k) => k[format])
    .map((k) => ({ name: downloadName(k, format), bytes: dataUrlBytes(k[format]) }));
  const zipName = TITLE.replace(/\.[^.]*$/, '') + '_karyograms_' + format + '.zip';
  saveFile(makeZip(files), zipName, 'application/zip');
};

window.addEventListener('resize', () => { if (fitMode) fit(); else setScale(scale); });

show(0);
</script>
</body>
</html>
)HTML";

    html.close();

    if (!html)
    {
        std::cerr << "ERROR: Failed to write web viewer file " << htmlPath << "\n";
        return false;
    }

    return true;
}


bool writeZipArchive(
    const std::vector<std::pair<std::string, std::vector<unsigned char>>>& entries,
    const std::filesystem::path& zipPath)
{
    // CRC-32 checksum the zip format stores for every file, so unzip tools can check the contents.
    const auto crc32 = [](const std::vector<unsigned char>& bytes) -> std::uint32_t
    {
        std::uint32_t crc = 0xFFFFFFFFu;

        for (unsigned char byte : bytes)
        {
            crc ^= byte;

            for (int bit = 0; bit < 8; ++bit)
            {
                crc = (crc >> 1) ^ (0xEDB88320u & (0u - (crc & 1u)));
            }
        }

        return crc ^ 0xFFFFFFFFu;
    };

    // Zip numbers are little-endian, written byte by byte so the file is the same on every machine.
    const auto put16 = [](std::vector<unsigned char>& out, std::uint32_t value)
    {
        out.push_back(static_cast<unsigned char>(value & 0xFF));
        out.push_back(static_cast<unsigned char>((value >> 8) & 0xFF));
    };
    const auto put32 = [&put16](std::vector<unsigned char>& out, std::uint32_t value)
    {
        put16(out, value & 0xFFFF);
        put16(out, value >> 16);
    };

    // File date and time in the zip's MS-DOS format, so the extracted files show when they were made.
    const std::time_t now = std::time(nullptr);
    const std::tm* local = std::localtime(&now);
    const std::uint32_t dosTime = local ? static_cast<std::uint32_t>((local->tm_hour << 11) | (local->tm_min << 5) | (local->tm_sec / 2)) : 0;
    const std::uint32_t dosDate = local ? static_cast<std::uint32_t>(((local->tm_year - 80) << 9) | ((local->tm_mon + 1) << 5) | local->tm_mday) : 0x21;

    std::vector<unsigned char> archive;
    std::vector<unsigned char> centralDirectory;

    for (const auto& [name, data] : entries)
    {
        // Safeguard: this simple writer has no ZIP64 support, so each file and the whole archive must stay under 4 GB.
        if (data.size() > 0xFFFFFFFFu || archive.size() > 0xFFFFFFFFu || name.size() > 0xFFFFu)
        {
            std::cerr << "ERROR: " << name << " is too large to add to the zip file " << zipPath << "\n";
            return false;
        }

        const std::uint32_t crc = crc32(data);
        const std::uint32_t size = static_cast<std::uint32_t>(data.size());
        const std::uint32_t offset = static_cast<std::uint32_t>(archive.size());

        put32(archive, 0x04034B50u);                                            // Local file header
        put16(archive, 20);                                                     // Version needed to extract (2.0)
        put16(archive, 0);                                                      // Flags
        put16(archive, 0);                                                      // Stored, no compression
        put16(archive, dosTime);
        put16(archive, dosDate);
        put32(archive, crc);
        put32(archive, size);                                                   // Compressed size (same, stored)
        put32(archive, size);
        put16(archive, static_cast<std::uint32_t>(name.size()));
        put16(archive, 0);                                                      // No extra field
        archive.insert(archive.end(), name.begin(), name.end());
        archive.insert(archive.end(), data.begin(), data.end());

        put32(centralDirectory, 0x02014B50u);                                   // Central directory entry for the same file
        put16(centralDirectory, 20);                                            // Made by version 2.0
        put16(centralDirectory, 20);
        put16(centralDirectory, 0);
        put16(centralDirectory, 0);
        put16(centralDirectory, dosTime);
        put16(centralDirectory, dosDate);
        put32(centralDirectory, crc);
        put32(centralDirectory, size);
        put32(centralDirectory, size);
        put16(centralDirectory, static_cast<std::uint32_t>(name.size()));
        put16(centralDirectory, 0);                                             // No extra field
        put16(centralDirectory, 0);                                             // No comment
        put16(centralDirectory, 0);                                             // Disk number
        put16(centralDirectory, 0);                                             // Internal attributes
        put32(centralDirectory, 0);                                             // External attributes
        put32(centralDirectory, offset);
        centralDirectory.insert(centralDirectory.end(), name.begin(), name.end());
    }

    if (entries.size() > 0xFFFFu)
    {
        std::cerr << "ERROR: Too many files for the zip file " << zipPath << "\n";
        return false;
    }

    const std::uint32_t directoryOffset = static_cast<std::uint32_t>(archive.size());
    archive.insert(archive.end(), centralDirectory.begin(), centralDirectory.end());

    put32(archive, 0x06054B50u);                                                // End of central directory record
    put16(archive, 0);
    put16(archive, 0);
    put16(archive, static_cast<std::uint32_t>(entries.size()));
    put16(archive, static_cast<std::uint32_t>(entries.size()));
    put32(archive, static_cast<std::uint32_t>(centralDirectory.size()));
    put32(archive, directoryOffset);
    put16(archive, 0);                                                          // No archive comment

    std::ofstream zipFile(zipPath, std::ios::binary);
    zipFile.write(reinterpret_cast<const char*>(archive.data()), static_cast<std::streamsize>(archive.size()));

    if (!zipFile)
    {
        std::cerr << "ERROR: Failed to write zip file " << zipPath << "\n";
        return false;
    }

    return true;
}


bool canOpenBrowser()
{
#if defined(_WIN32)
    return true;                                                                // Native Windows always has a desktop and a default browser.
#elif defined(__APPLE__)
    return !isRemoteSSHsession();                                               // A Mac reached over SSH would open the browser on its own (unseen) screen.
#else
    if (isRunningInWSL())
    {
        return wslCanLaunchWindowsPrograms();                                   // WSL opens the Windows browser, which needs Windows interop.
    }

    // Plain Linux: a browser can only appear if there is a graphical display. Desktops set DISPLAY (X11) or
    // WAYLAND_DISPLAY (Wayland); SSH with X forwarding (ssh -X) also sets DISPLAY, so that case still works.
    // Why: without a display, xdg-open either fails or, worse, starts a text-mode browser (such as lynx or w3m)
    // inside the terminal, which takes over the terminal and blocks the program until the user quits it.
    return std::getenv("DISPLAY") != nullptr || std::getenv("WAYLAND_DISPLAY") != nullptr;
#endif
}


bool openInBrowser(const std::filesystem::path& htmlPath)
{
    if (!canOpenBrowser())                                                      // Safeguard for headless machines, see canOpenBrowser().
    {
        return false;
    }

    const std::string path = std::filesystem::absolute(htmlPath).string();

#if defined(_WIN32)
    // Windows file names cannot contain double quotes, so a quote here means something is wrong with the path.
    // Refusing it keeps the quoting of the "start" command below intact.
    if (path.find('"') != std::string::npos)
    {
        return false;
    }
    const std::string command = "start \"\" \"" + path + "\"";                    // The empty "" is the window title, start would otherwise treat the quoted path as the title.
#elif defined(__APPLE__)
    const std::string command = "open " + shellQuote(path);
#else
    if (isRunningInWSL())
    {
        // wslpath -w turns /tmp/x.html into \\wsl.localhost\<distro>\tmp\x.html, a path Windows programs can read,
        // and explorer.exe opens it with the user's default Windows browser.
        const std::string wslCommand =
            "explorer.exe \"$(wslpath -w " + shellQuote(path) + ")\" > /dev/null 2>&1";

        // explorer.exe returns 1 even when it succeeds, so its exit code cannot tell us anything. It is stored in a
        // variable only because some compilers warn when the result of std::system() is ignored.
        const int explorerStatus = std::system(wslCommand.c_str());
        (void)explorerStatus;
        return true;
    }

    const std::string command = "xdg-open " + shellQuote(path) + " > /dev/null 2>&1";   // Output is hidden because xdg-open can print browser noise into the CLI output.
#endif

    return std::system(command.c_str()) == 0;
}

} // namespace sddparser
