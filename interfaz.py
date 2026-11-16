import tkinter as tk
from tkinter import ttk, filedialog, messagebox


def hash_function(key, a, length):
    h = 0
    for ch in key:
        h = (h * a + ord(ch)) % length
    return h


def es_primo(n):
    if n < 2:
        return False
    i = 2
    while i * i <= n:
        if n % i == 0:
            return False
        i += 1
    return True


def siguiente_primo(n):
    while not es_primo(n):
        n += 1
    return n


def leer_fasta(ruta):
    """Devuelve una lista de [id, secuencia]."""
    lista, actual = [], None
    with open(ruta, encoding="utf-8", errors="replace") as f:
        for linea in f:
            linea = linea.strip()
            if not linea:
                continue
            if linea[0] == ">":
                if actual:
                    lista.append(actual)
                actual = [linea[1:], ""]
            elif actual is not None:
                actual[1] += linea
    if actual:
        lista.append(actual)
    return lista


class TablaHash:
    """Tabla 1: linear probing, se agranda cuando el factor de carga pasa de 0.7"""

    def __init__(self, tam):
        self.crear(tam)

    def crear(self, tam):
        self.length = tam
        self.num_keys = 0
        self.bins = [None] * tam          # cada bin: None o [id, secuencia]

    def busqueda(self, key):
        index = hash_function(key, 37, self.length)
        conteo = 0
        while self.bins[index] is not None and conteo != self.length:
            if self.bins[index][0] == key:
                return index
            index = (index + 1) % self.length
            conteo += 1
        return -1

    def redimensionar(self):
        vieja = self.bins
        self.crear(siguiente_primo(2 * (self.num_keys + 1)))
        for b in vieja:
            if b is not None:
                i = hash_function(b[0], 37, self.length)
                while self.bins[i] is not None:
                    i = (i + 1) % self.length
                self.bins[i] = b
                self.num_keys += 1

    def insertar(self, id_, seq):
        pos = self.busqueda(id_)
        if pos != -1:
            self.bins[pos] = [id_, seq]
            return True
        if (self.num_keys + 1) / self.length > 0.7:
            self.redimensionar()
        index = hash_function(id_, 37, self.length)
        conteo = 0
        while self.bins[index] is not None and conteo != self.length:
            index = (index + 1) % self.length
            conteo += 1
        if conteo == self.length:
            return False
        self.bins[index] = [id_, seq]
        self.num_keys += 1
        return True


class TablaKmer:
    """Tabla 2: double hashing, 4^k bins"""

    def __init__(self, k):
        self.k = k
        self.length = 4 ** k
        self.num_keys = 0
        self.bins = [None] * self.length  # cada bin: None o [kmer, veces]

    def _h1(self, key):
        return hash_function(key, 37, self.length)

    def _h2(self, key):
        return 2 * hash_function(key, 41, self.length // 2) + 1

    def insertar(self, key):
        index, salto, conteo = self._h1(key), self._h2(key), 0
        while conteo != self.length:
            if self.bins[index] is None:
                self.bins[index] = [key, 1]
                self.num_keys += 1
                return True
            if self.bins[index][0] == key:
                self.bins[index][1] += 1
                return True
            index = (index + salto) % self.length
            conteo += 1
        return False

    def busqueda(self, key):
        index, salto, conteo = self._h1(key), self._h2(key), 0
        while self.bins[index] is not None and conteo != self.length:
            if self.bins[index][0] == key:
                return self.bins[index][1]
            index = (index + salto) % self.length
            conteo += 1
        return 0

    def construir(self, seq):
        k = self.k
        for i in range(len(seq) - k + 1):
            kmer = seq[i:i + k]
            if all(c in "ACGT" for c in kmer):
                self.insertar(kmer)


# ---------- Max-Heap ----------
def heapify_down(heap, tamano, i):
    while True:
        mayor, L, R = i, 2 * i + 1, 2 * i + 2
        if L < tamano and heap[L][1] > heap[mayor][1]:
            mayor = L
        if R < tamano and heap[R][1] > heap[mayor][1]:
            mayor = R
        if mayor == i:
            return
        heap[i], heap[mayor] = heap[mayor], heap[i]
        i = mayor


def construir_heap(heap):
    for i in range((len(heap) - 1) // 2, -1, -1):
        heapify_down(heap, len(heap), i)


def extraer_max(heap):
    mayor = heap[0]
    heap[0] = heap[-1]
    heap.pop()
    if heap:
        heapify_down(heap, len(heap), 0)
    return mayor


def top_n_kmers(tabla, n):
    heap = [[b[0], b[1]] for b in tabla.bins if b is not None]
    construir_heap(heap)
    salida = []
    while len(salida) < n and heap:
        salida.append(extraer_max(heap))
    return salida


# ---------- Matriz de adyacencia ----------
def matriz_adyacencia(secuencia, beta):
    idx = {"A": 0, "C": 1, "G": 2, "T": 3}
    m = len(secuencia) // 2
    par, usado = [], set()
    for p in range(m):
        a, b = idx.get(secuencia[2 * p], -1), idx.get(secuencia[2 * p + 1], -1)
        if a == -1 or b == -1:
            par.append(-1)
        else:
            par.append(a * 4 + b)
            usado.add(a * 4 + b)
    A = [[0.0] * 16 for _ in range(16)]
    for i in range(m):
        if par[i] == -1:
            continue
        for j in range(i + 1, m):
            if par[j] != -1:
                A[par[i]][par[j]] += 1.0 / ((j - i) ** beta)
    return sorted(usado), A


def nombre_2mer(c):
    return "ACGT"[c // 4] + "ACGT"[c % 4]


# =====================================================================
#  INTERFAZ
# =====================================================================
BG, PANEL, CARD = "#0f172a", "#1e293b", "#334155"
TXT, MUTED = "#e2e8f0", "#94a3b8"
ACC, ACC2, OK, WARN = "#8b5cf6", "#22d3ee", "#4ade80", "#fbbf24"
BASES = {"A": "#4ade80", "C": "#60a5fa", "G": "#fbbf24", "T": "#f87171"}


def lerp(c1, c2, t):
    a = [int(c1[i:i + 2], 16) for i in (1, 3, 5)]
    b = [int(c2[i:i + 2], 16) for i in (1, 3, 5)]
    return "#%02x%02x%02x" % tuple(int(a[i] + (b[i] - a[i]) * t) for i in range(3))


def boton(parent, texto, cmd, color=ACC, fg="white"):
    return tk.Button(parent, text=texto, command=cmd, bg=color, fg=fg,
                     activebackground=lerp(color, "#ffffff", 0.2), activeforeground=fg,
                     relief="flat", bd=0, font=("Segoe UI", 10, "bold"),
                     padx=14, pady=6, cursor="hand2")


def etiqueta(parent, texto, fg=TXT, bg=BG, size=10, bold=False):
    return tk.Label(parent, text=texto, bg=bg, fg=fg,
                    font=("Segoe UI", size, "bold" if bold else "normal"))


def entrada(parent, ancho=20):
    return tk.Entry(parent, bg=CARD, fg=TXT, insertbackground=TXT, relief="flat",
                    font=("Segoe UI", 11), width=ancho)


class App:
    def __init__(self, root):
        self.root = root
        root.title("DNA Lab")
        root.geometry("960x760")
        root.configure(bg=BG)
        self.ht = TablaHash(10)
        self.ids = []
        self.hk = None
        self.id_actual = None
        self.seq = ""
        self.estilos()
        self.cabecera()
        self.nb = ttk.Notebook(root)
        self.nb.pack(fill="both", expand=True, padx=18, pady=(6, 16))
        self.tab_secuencias()
        self.tab_kmers()
        self.tab_top()
        self.tab_matriz()

    # ---------- estilos y piezas ----------
    def estilos(self):
        s = ttk.Style()
        s.theme_use("clam")
        s.configure("TNotebook", background=BG, borderwidth=0)
        s.configure("TNotebook.Tab", background=PANEL, foreground=MUTED,
                    padding=(18, 8), font=("Segoe UI", 11, "bold"), borderwidth=0)
        s.map("TNotebook.Tab", background=[("selected", ACC)],
              foreground=[("selected", "white")])
        s.configure("TCombobox", fieldbackground=CARD, background=CARD,
                    foreground=TXT, arrowcolor=TXT, borderwidth=0)
        s.map("TCombobox", fieldbackground=[("readonly", CARD)],
              foreground=[("readonly", TXT)])
        s.configure("Treeview", background=PANEL, fieldbackground=PANEL,
                    foreground=TXT, rowheight=24, borderwidth=0)
        s.configure("Treeview.Heading", background=CARD, foreground=TXT,
                    font=("Segoe UI", 10, "bold"), borderwidth=0)
        s.map("Treeview", background=[("selected", ACC)])
        self.root.option_add("*TCombobox*Listbox.background", PANEL)
        self.root.option_add("*TCombobox*Listbox.foreground", TXT)

    def tarjeta(self, parent, titulo, color):
        f = tk.Frame(parent, bg=PANEL)
        f.pack(side="left", padx=(0, 12))
        tk.Frame(f, bg=color, width=5).pack(side="left", fill="y")
        interior = tk.Frame(f, bg=PANEL)
        interior.pack(side="left", padx=14, pady=8)
        etiqueta(interior, titulo, MUTED, PANEL, 9).pack(anchor="w")
        valor = etiqueta(interior, "-", color, PANEL, 20, True)
        valor.pack(anchor="w")
        return valor

    def nuevo_tab(self, titulo):
        f = tk.Frame(self.nb, bg=BG, padx=14, pady=8)
        self.nb.add(f, text=titulo)
        return f

    def cabecera(self):
        top = tk.Frame(self.root, bg=BG)
        top.pack(fill="x", padx=18, pady=(14, 6))
        etiqueta(top, "🧬 DNA Lab", ACC2, BG, 22, True).pack(side="left")
        etiqueta(top, "   hash · k-mers · heap · matriz", MUTED, BG, 10).pack(side="left", pady=(10, 0))
        boton(top, "📂  Abrir archivo FASTA", self.abrir).pack(side="right")

        barra = tk.Frame(self.root, bg=PANEL)
        barra.pack(fill="x", padx=18, pady=6)
        etiqueta(barra, "Secuencia activa:", TXT, PANEL, 10, True).pack(side="left", padx=(12, 8), pady=10)
        self.combo = ttk.Combobox(barra, state="readonly", width=48)
        self.combo.pack(side="left")
        self.combo.bind("<<ComboboxSelected>>", lambda e: self.elegir(self.combo.get()))
        self.estado = etiqueta(barra, "👈 Empieza abriendo un archivo FASTA", MUTED, PANEL)
        self.estado.pack(side="right", padx=12)

    def requiere_seq(self):
        if not self.seq:
            messagebox.showinfo("Falta algo", "Primero abre un FASTA y elige una secuencia.")
            return False
        return True

    def requiere_kmers(self):
        if self.hk is None:
            messagebox.showinfo("Falta algo", "Primero construye la tabla en la pestaña 'K-mers'.")
            return False
        return True

    # ---------- pestaña 1: secuencias ----------
    def tab_secuencias(self):
        f = self.nuevo_tab("🔎 Secuencias")
        fila = tk.Frame(f, bg=BG)
        fila.pack(fill="x", pady=(10, 6))
        etiqueta(fila, "Buscar ID en la tabla hash:").pack(side="left")
        self.e_id = entrada(fila, 28)
        self.e_id.pack(side="left", padx=8, ipady=4)
        boton(fila, "Buscar", self.buscar_id, ACC2, "#0f172a").pack(side="left")
        self.r_busq = etiqueta(fila, "", OK)
        self.r_busq.pack(side="left", padx=12)

        cards = tk.Frame(f, bg=BG)
        cards.pack(fill="x", pady=10)
        self.c_keys = self.tarjeta(cards, "Secuencias guardadas", OK)
        self.c_bins = self.tarjeta(cards, "Bins de la tabla", ACC2)
        self.c_load = self.tarjeta(cards, "Factor de carga", WARN)
        self.c_len = self.tarjeta(cards, "Largo de la secuencia", ACC)

        etiqueta(f, "Composición de bases", MUTED, BG, 9).pack(anchor="w")
        self.c_comp = tk.Canvas(f, height=26, bg=BG, highlightthickness=0)
        self.c_comp.pack(fill="x", pady=(2, 10))

        self.t_seq = tk.Text(f, bg=PANEL, fg=TXT, font=("Consolas", 12), relief="flat",
                             wrap="char", padx=12, pady=12, state="disabled")
        self.t_seq.pack(fill="both", expand=True)
        for b, c in BASES.items():
            self.t_seq.tag_configure(b, foreground=c)
        self.t_seq.tag_configure("titulo", foreground=MUTED)

    def abrir(self):
        ruta = filedialog.askopenfilename(
            filetypes=[("FASTA", "*.fasta *.fa *.fna *.txt"), ("Todos", "*.*")])
        if not ruta:
            return
        try:
            lista = leer_fasta(ruta)
        except OSError:
            return messagebox.showerror("Error", "No se pudo abrir el archivo.")
        if not lista:
            return messagebox.showwarning("Aviso", "No encontré secuencias (¿empiezan con '>'?).")
        for id_, seq in lista:
            self.ht.insertar(id_, seq)
            if id_ not in self.ids:
                self.ids.append(id_)
        self.combo["values"] = self.ids
        self.estado.config(text=f"✅ {len(lista)} secuencias leídas", fg=OK)
        self.elegir(lista[0][0])

    def elegir(self, id_):
        pos = self.ht.busqueda(id_)            # se busca de verdad en la tabla hash
        if pos == -1:
            return
        self.id_actual, self.seq = self.ht.bins[pos]
        self.combo.set(id_)
        self.hk = None
        self.limpiar_kmers()
        self.actualizar_secuencia()

    def buscar_id(self):
        id_ = self.e_id.get().strip()
        pos = self.ht.busqueda(id_)
        if pos == -1:
            self.r_busq.config(text="❌ No encontrada", fg="#f87171")
        else:
            self.r_busq.config(text=f"✅ Encontrada en el bin {pos}", fg=OK)
            self.elegir(id_)

    def actualizar_secuencia(self):
        self.c_keys.config(text=str(self.ht.num_keys))
        self.c_bins.config(text=str(self.ht.length))
        self.c_load.config(text=f"{self.ht.num_keys / self.ht.length:.2f}")
        self.c_len.config(text=f"{len(self.seq):,}")

        t = self.t_seq
        t.config(state="normal")
        t.delete("1.0", "end")
        t.insert("end", f">{self.id_actual}\n\n", "titulo")
        for ch in self.seq[:1500]:
            t.insert("end", ch, ch if ch in BASES else ())
        if len(self.seq) > 1500:
            t.insert("end", "\n\n… (se muestran las primeras 1500 letras)", "titulo")
        t.config(state="disabled")

        c = self.c_comp
        c.delete("all")
        c.update_idletasks()
        ancho = max(c.winfo_width(), 600)
        total = sum(self.seq.count(b) for b in BASES)
        x = 0
        for b, color in BASES.items():
            if total == 0:
                break
            frac = self.seq.count(b) / total
            w = frac * ancho
            c.create_rectangle(x, 0, x + w, 26, fill=color, outline="")
            if w > 70:
                c.create_text(x + w / 2, 13, text=f"{b}  {frac * 100:.1f}%",
                              fill="#0f172a", font=("Segoe UI", 10, "bold"))
            x += w

    # ---------- pestaña 2: k-mers ----------
    def tab_kmers(self):
        f = self.nuevo_tab("🔤 K-mers")
        fila = tk.Frame(f, bg=BG)
        fila.pack(fill="x", pady=(10, 6))
        etiqueta(fila, "Tamaño k:").pack(side="left")
        self.s_k = tk.Scale(fila, from_=1, to=10, orient="horizontal", length=220,
                            bg=BG, fg=TXT, troughcolor=CARD, highlightthickness=0,
                            activebackground=ACC, font=("Segoe UI", 10, "bold"))
        self.s_k.set(3)
        self.s_k.pack(side="left", padx=10)
        boton(fila, "⚙ Construir tabla", self.construir).pack(side="left")

        cards = tk.Frame(f, bg=BG)
        cards.pack(fill="x", pady=10)
        self.c_kbins = self.tarjeta(cards, "Bins (4^k)", ACC2)
        self.c_kdist = self.tarjeta(cards, "K-mers distintos", OK)

        fila2 = tk.Frame(f, bg=BG)
        fila2.pack(fill="x", pady=6)
        etiqueta(fila2, "¿Cuántas veces aparece?").pack(side="left")
        self.e_kmer = entrada(fila2, 14)
        self.e_kmer.pack(side="left", padx=8, ipady=4)
        boton(fila2, "Buscar k-mer", self.buscar_kmer, ACC2, "#0f172a").pack(side="left")
        self.r_kmer = etiqueta(fila2, "", OK, BG, 11, True)
        self.r_kmer.pack(side="left", padx=12)

        etiqueta(f, "Primeros 20 k-mers (en el orden de la tabla)", MUTED, BG, 9).pack(anchor="w", pady=(10, 2))
        self.tree = ttk.Treeview(f, columns=("bin", "kmer", "veces"), show="headings", height=12)
        for col, txt, w in (("bin", "Bin", 90), ("kmer", "K-mer", 200), ("veces", "Veces", 120)):
            self.tree.heading(col, text=txt)
            self.tree.column(col, width=w, anchor="center")
        self.tree.pack(fill="both", expand=True)

    def limpiar_kmers(self):
        self.c_kbins.config(text="-")
        self.c_kdist.config(text="-")
        self.r_kmer.config(text="")
        self.tree.delete(*self.tree.get_children())

    def construir(self):
        if not self.requiere_seq():
            return
        self.root.config(cursor="watch")
        self.root.update()
        self.hk = TablaKmer(self.s_k.get())
        self.hk.construir(self.seq)
        self.root.config(cursor="")
        self.c_kbins.config(text=f"{self.hk.length:,}")
        self.c_kdist.config(text=f"{self.hk.num_keys:,}")
        self.tree.delete(*self.tree.get_children())
        mostrados = 0
        for i, b in enumerate(self.hk.bins):
            if b is not None:
                self.tree.insert("", "end", values=(i, b[0], b[1]))
                mostrados += 1
                if mostrados == 20:
                    break

    def buscar_kmer(self):
        if not self.requiere_kmers():
            return
        km = self.e_kmer.get().strip().upper()
        if len(km) != self.hk.k or any(c not in "ACGT" for c in km):
            return self.r_kmer.config(text=f"Escribe {self.hk.k} letras (A, C, G, T)", fg="#f87171")
        self.r_kmer.config(text=f"{km} aparece {self.hk.busqueda(km)} veces", fg=OK)

    # ---------- pestaña 3: top n (max-heap) ----------
    def tab_top(self):
        f = self.nuevo_tab("🏆 Top k-mers")
        fila = tk.Frame(f, bg=BG)
        fila.pack(fill="x", pady=(10, 6))
        etiqueta(fila, "¿Cuántos quieres ver?").pack(side="left")
        self.sp_n = tk.Spinbox(fila, from_=1, to=50, width=5, bg=CARD, fg=TXT,
                               buttonbackground=CARD, relief="flat",
                               font=("Segoe UI", 12), insertbackground=TXT)
        self.sp_n.delete(0, "end")
        self.sp_n.insert(0, "5")
        self.sp_n.pack(side="left", padx=8, ipady=3)
        boton(fila, "🏆 Ver top", self.ver_top).pack(side="left")

        cont = tk.Frame(f, bg=PANEL)
        cont.pack(fill="both", expand=True, pady=10)
        self.cv_top = tk.Canvas(cont, bg=PANEL, highlightthickness=0)
        sb = tk.Scrollbar(cont, command=self.cv_top.yview)
        self.cv_top.config(yscrollcommand=sb.set)
        sb.pack(side="right", fill="y")
        self.cv_top.pack(side="left", fill="both", expand=True)

    def ver_top(self):
        if not self.requiere_kmers():
            return
        try:
            n = int(self.sp_n.get())
            assert n >= 1
        except (ValueError, AssertionError):
            return messagebox.showwarning("Aviso", "n debe ser un entero mayor que 0.")
        datos = top_n_kmers(self.hk, n)
        c = self.cv_top
        c.delete("all")
        if not datos:
            return
        mx = datos[0][1]
        for i, (km, cn) in enumerate(datos):
            y = 14 + i * 36
            c.create_text(14, y + 14, text=f"{i + 1}.", anchor="w", fill=MUTED, font=("Segoe UI", 11, "bold"))
            for j, ch in enumerate(km):
                c.create_text(52 + j * 12, y + 14, text=ch, anchor="w",
                              fill=BASES.get(ch, TXT), font=("Consolas", 14, "bold"))
            ancho = max(int(520 * cn / mx), 4)
            color = lerp(ACC2, ACC, i / max(1, len(datos) - 1))
            c.create_rectangle(190, y, 190 + ancho, y + 28, fill=color, outline="")
            c.create_text(190 + ancho + 10, y + 14, text=f"{cn:,}", anchor="w",
                          fill=TXT, font=("Segoe UI", 11, "bold"))
        c.config(scrollregion=(0, 0, 900, 14 + len(datos) * 36 + 10))

    # ---------- pestaña 4: matriz de adyacencia ----------
    def tab_matriz(self):
        f = self.nuevo_tab("🧩 Matriz")
        fila = tk.Frame(f, bg=BG)
        fila.pack(fill="x", pady=(10, 6))
        etiqueta(fila, "beta:").pack(side="left")
        self.e_beta = entrada(fila, 6)
        self.e_beta.insert(0, "1")
        self.e_beta.pack(side="left", padx=8, ipady=4)
        boton(fila, "🧩 Calcular matriz", self.calcular_matriz).pack(side="left")
        self.l_hover = etiqueta(fila, "Pasa el mouse por una celda 👆", MUTED)
        self.l_hover.pack(side="left", padx=16)
        self.cv_mat = tk.Canvas(f, bg=PANEL, highlightthickness=0)
        self.cv_mat.pack(fill="both", expand=True, pady=10)

    def calcular_matriz(self):
        if not self.requiere_seq():
            return
        try:
            beta = float(self.e_beta.get().replace(",", "."))
        except ValueError:
            return messagebox.showwarning("Aviso", "beta debe ser un número.")
        usados, A = matriz_adyacencia(self.seq, beta)
        c = self.cv_mat
        c.delete("all")
        n = len(usados)
        if n == 0:
            return c.create_text(30, 30, text="La secuencia es muy corta o no tiene A, C, G, T.",
                                 anchor="w", fill=MUTED, font=("Segoe UI", 12))
        mx = max(A[a][b] for a in usados for b in usados) or 1
        cell = max(30, min(52, 600 // (n + 1)))
        ox = oy = cell
        fuente = ("Consolas", 7 if cell < 40 else 8)
        for k, u in enumerate(usados):
            nom = nombre_2mer(u)
            for pos, (x, y) in enumerate(((ox + k * cell + cell / 2, oy / 2),
                                          (ox / 2, oy + k * cell + cell / 2))):
                c.create_text(x, y, text=nom, fill=TXT, font=("Consolas", 9, "bold"))
        for i, a in enumerate(usados):
            for j, b in enumerate(usados):
                v = A[a][b]
                t = v / mx
                color = lerp("#1e1b4b", "#f472b6", t ** 0.6) if v else "#172033"
                x0, y0 = ox + j * cell, oy + i * cell
                r = c.create_rectangle(x0 + 1, y0 + 1, x0 + cell - 1, y0 + cell - 1,
                                       fill=color, outline="")
                if v:
                    c.create_text(x0 + cell / 2, y0 + cell / 2, text=f"{v:.2f}",
                                  fill="white" if t > 0.35 else MUTED, font=fuente)
                txt = f"{nombre_2mer(a)} → {nombre_2mer(b)}:  {v:.3f}"
                c.tag_bind(r, "<Enter>", lambda e, s=txt: self.l_hover.config(text=s, fg=ACC2))


if __name__ == "__main__":
    try:  # pantallas de alta resolucion en Windows (texto mas nitido)
        import ctypes
        ctypes.windll.shcore.SetProcessDpiAwareness(1)
    except Exception:
        pass
    root = tk.Tk()
    App(root)
    root.mainloop()