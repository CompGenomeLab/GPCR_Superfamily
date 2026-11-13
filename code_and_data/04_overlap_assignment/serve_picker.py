"""Serve overlap_mrca_picker.html and auto-save its selections.

Run:  python serve_picker.py   (after build_picker.py has written overlap_mrca_picker.html here)

The page POSTs to /save on every click, and this writes mrca_selections.tsv
next to it, so nothing has to be exported or clicked by hand.
"""
import http.server, socketserver, os, sys, webbrowser, threading, datetime

HERE = os.path.dirname(os.path.abspath(__file__))
OUT  = os.path.join(HERE, "mrca_selections.tsv")
PORT = 8777


class Handler(http.server.SimpleHTTPRequestHandler):
    def __init__(self, *a, **kw):
        super().__init__(*a, directory=HERE, **kw)

    def do_POST(self):
        if self.path.rstrip("/") != "/save":
            self.send_error(404)
            return
        n = int(self.headers.get("Content-Length", 0))
        body = self.rfile.read(n).decode("utf-8")
        # keep one backup of the previous state before overwriting
        if os.path.exists(OUT):
            try:
                with open(OUT, encoding="utf-8") as fh:
                    prev = fh.read()
                if prev != body:
                    with open(OUT + ".bak", "w", encoding="utf-8", newline="\n") as fh:
                        fh.write(prev)
            except OSError:
                pass
        with open(OUT, "w", encoding="utf-8", newline="\n") as fh:
            fh.write(body)
        rows = max(body.count("\n") - 1, 0)
        print("  %s  saved %d rows -> mrca_selections.tsv"
              % (datetime.datetime.now().strftime("%H:%M:%S"), rows))
        self.send_response(204)
        self.send_header("Access-Control-Allow-Origin", "*")
        self.end_headers()

    def end_headers(self):
        self.send_header("Cache-Control", "no-store")
        super().end_headers()

    def log_message(self, *a):
        pass


if __name__ == "__main__":
    url = "http://localhost:%d/overlap_mrca_picker.html" % PORT
    socketserver.TCPServer.allow_reuse_address = True
    try:
        srv = socketserver.TCPServer(("127.0.0.1", PORT), Handler)
    except OSError as e:
        print("Could not start on port %d (%s)." % (PORT, e))
        print("It may already be running -- just open %s" % url)
        sys.exit(1)
    print("MRCA picker running.")
    print("  %s" % url)
    print("  saving to %s" % OUT)
    print("  leave this window open; Ctrl+C to stop\n")
    threading.Timer(0.6, lambda: webbrowser.open(url)).start()
    try:
        srv.serve_forever()
    except KeyboardInterrupt:
        print("\nstopped")
