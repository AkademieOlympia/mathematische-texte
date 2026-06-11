import math

try:
    import tkinter as tk
    from tkinter import messagebox, simpledialog
except Exception:
    tk = None
    messagebox = None
    simpledialog = None


def ask_number():
    """Fragt eine Zahl N per Dialog ab, mit Konsolen-Fallback."""
    if tk is not None and simpledialog is not None:
        root = tk.Tk()
        root.withdraw()
        try:
            value = simpledialog.askfloat("Eingabe", "Bitte Zahl N eingeben:")
        finally:
            root.destroy()
        if value is None:
            raise SystemExit("Abgebrochen.")
        return value

    raw = input("Bitte Zahl N eingeben: ").strip()
    if not raw:
        raise SystemExit("Keine Zahl eingegeben.")
    return float(raw)


def main():
    n = ask_number()
    if n < 0:
        raise SystemExit("N muss >= 0 sein.")

    fourth_root = math.sqrt(math.sqrt(n))
    integers = list(range(int(fourth_root), 0, -1))

    print(f"N = {n}")
    print(f"Vierte Wurzel von N: {fourth_root}")
    print("Ganzzahlen abwaerts von int(vierte Wurzel) bis 1:")
    for value in integers:
        print(value)

    if tk is not None and messagebox is not None:
        root = tk.Tk()
        root.withdraw()
        try:
            integer_text = ", ".join(str(value) for value in integers) if integers else "keine"
            messagebox.showinfo(
                "Ergebnis",
                f"N = {n}\n"
                f"Vierte Wurzel: {fourth_root}\n"
                f"Ganzzahlen abwaerts bis 1: {integer_text}",
            )
        finally:
            root.destroy()


if __name__ == "__main__":
    main()
