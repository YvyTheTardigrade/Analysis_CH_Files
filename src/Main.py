import GUI  
import tkinter as tk
import multiprocessing as mp



def main() -> int:
    mp.freeze_support()
    root = tk.Tk()
    app = GUI.App(root)
    root.mainloop()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
