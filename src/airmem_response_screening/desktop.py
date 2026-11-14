"""Single-purpose desktop front end; no network service or browser."""
from datetime import datetime
from pathlib import Path
import json
from . import __version__
from queue import Empty, Queue
from threading import Thread
import tkinter as tk
from tkinter import filedialog, messagebox, ttk

from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
from matplotlib.figure import Figure

from .analysis import resolve_monitoring, run
from .figures import TITLES, build, export, read_run


class Window:
    def __init__(self, root):
        self.root=root; root.title(f'AirMem Response-Side Screening v{__version__}'); root.geometry('1280x780')
        root.minsize(930,610); root.protocol('WM_DELETE_WINDOW',self.close)
        self.data=self.original=self.output=self.job=None
        self.queue=Queue()
        frame=ttk.Frame(root,padding=12); frame.pack(fill='both',expand=True)
        actions=ttk.Frame(frame);actions.pack(fill='x',pady=(0,8))
        self.select=ttk.Button(actions,text='Select monitoring data…',command=self.choose_data)
        self.choose_time=ttk.Button(actions,text='Optional time-vector ZIP…',command=self.choose_original)
        self.calculate=ttk.Button(actions,text='Run monitoring calculation…',command=self.calculate_run)
        self.resume=ttk.Button(actions,text='Open completed run…',command=self.open_run)
        for button in (self.select,self.choose_time,self.calculate,self.resume):
            button.pack(side='left',padx=(0,8))
        self.data_label=tk.StringVar(value='Select 2022Q4 or the parent folder containing 2022Q4.')
        self.time_label=tk.StringVar(value='Public calendar index (if available).')
        self.run_label=tk.StringVar(value='No calculation loaded.')
        for key,var in [('Data',self.data_label),('Figure 3 time source',self.time_label),('Run',self.run_label)]:
            row=ttk.Frame(frame);row.pack(fill='x',pady=1)
            ttk.Label(row,text=key+':',width=22).pack(side='left',anchor='n');ttk.Label(row,textvariable=var,wraplength=790).pack(side='left',fill='x',expand=True,anchor='w')
        pane=ttk.Panedwindow(frame,orient='horizontal');pane.pack(fill='both',expand=True,pady=10)
        left=ttk.Frame(pane,padding=(0,0,8,0));pane.add(left,weight=0)
        ttk.Label(left,text='V19.0 figures',padding=(4,4)).pack(anchor='w')
        items=ttk.Frame(left);items.pack(fill='both',expand=True)
        self.items=tk.Listbox(items,width=44,exportselection=False,activestyle='none')
        scroll=ttk.Scrollbar(items,orient='vertical',command=self.items.yview)
        self.items.configure(yscrollcommand=scroll.set);scroll.pack(side='right',fill='y');self.items.pack(side='left',fill='both',expand=True)
        self.numbers=list(TITLES)
        for n in self.numbers:self.items.insert('end',f'{n:02d}  {TITLES[n]}')
        self.items.selection_set(0);self.items.bind('<Double-1>',lambda _:self.draw_selected())
        self.draw=ttk.Button(left,text='Generate selected figure',command=self.draw_selected);self.draw.pack(fill='x',pady=(8,0))
        self.viewer=ttk.Frame(pane);pane.add(self.viewer,weight=1)
        self.message=tk.StringVar(value='2022Q4 analysis · English figures')
        ttk.Label(frame,textvariable=self.message,wraplength=1000).pack(fill='x',pady=(0,4))
        self.progress=ttk.Progressbar(frame,maximum=100);self.progress.pack(fill='x')
        self.canvas=self.navigation=None
        blank=Figure();ax=blank.subplots();ax.axis('off');ax.text(.5,.5,'Run the calculation, then choose a figure.',ha='center',va='center')
        self.show_figure(blank);self.refresh()

    def refresh(self):
        busy=self.job is not None
        for button in (self.select,self.choose_time,self.resume):button.configure(state='disabled' if busy else 'normal')
        self.calculate.configure(state='normal' if not busy and self.data else 'disabled')
        self.draw.configure(state='normal' if not busy and self.output else 'disabled')

    def choose_data(self):
        path=filedialog.askdirectory(parent=self.root,title='Select monitoring data')
        if not path:return
        try:self.data=resolve_monitoring(path)
        except Exception as exc:self.fail(str(exc));return
        self.data_label.set(str(self.data));self.refresh()

    def choose_original(self):
        if self.original is not None:
            self.original=None;self.time_label.set('Public calendar index (if available).')
            self.choose_time.configure(text='Optional time-vector ZIP…');return
        path=filedialog.askopenfilename(parent=self.root,title='Select original timestamp vectors',filetypes=[('Time-vector ZIP','*.zip')])
        if path:
            self.original=Path(path);self.time_label.set(str(self.original));self.choose_time.configure(text='Clear time-vector ZIP')

    def set_run(self,path):
        summary=json.loads((path/'summary.json').read_text(encoding='utf-8'))
        self.output=path
        self.run_label.set(f'{path}  |  {summary["window_count"]:,} windows; {summary["selected"]["MedianPeakFreq_Hz"]:.3f} Hz; {summary["paired"]:,} paired passes')
        self.message.set(f'Calculation ready. Figure 3: {summary["coverage_source"] or "unavailable"}.')

    def open_run(self):
        path=filedialog.askdirectory(parent=self.root,title='Select completed V19.0 run')
        if not path:return
        path=Path(path)
        try:
            read_run(path)
            self.set_run(path)
        except Exception as exc:self.fail(str(exc))
        self.refresh()

    def calculate_run(self):
        parent=filedialog.askdirectory(parent=self.root,title='Select parent folder for the new results')
        if not parent:return
        dest=Path(parent)/datetime.now().strftime('Response-V19-%Y%m%d-%H%M%S-%f')
        data,original=self.data,self.original
        self.run_job(lambda progress:run(data,dest,original_root=original,progress=progress),self.set_run)

    def draw_selected(self):
        if self.output is None or self.job is not None or not self.items.curselection():return
        n=self.numbers[self.items.curselection()[0]];output,data=self.output,self.data
        def draw(progress):
            progress(20,f'Calculating Figure {n:02d}');fig=build(n,output,data)
            progress(75,'Saving PNG and PDF');png,_=export(fig,n,output)
            return fig,png
        def display(result):
            self.show_figure(result[0]);self.message.set(f'Saved {result[1]}')
        self.run_job(draw,display)

    def show_figure(self,fig):
        for widget in (self.navigation,self.canvas.get_tk_widget() if self.canvas else None):
            if widget is not None:widget.destroy()
        self.canvas=FigureCanvasTkAgg(fig,master=self.viewer)
        self.navigation=NavigationToolbar2Tk(self.canvas,self.viewer,pack_toolbar=False)
        self.navigation.pack(side='top',fill='x');self.canvas.get_tk_widget().pack(fill='both',expand=True);self.canvas.draw_idle()

    def run_job(self,function,finished):
        if self.job is not None:return
        self.progress.configure(value=0)
        def work():
            try:self.queue.put(('done',function(lambda *info:self.queue.put(('progress',info)))))
            except Exception as exc:self.queue.put(('error',str(exc)))
        self.job=Thread(target=work,daemon=True);self.refresh();self.job.start()
        self.root.after(60,self.poll,finished)

    def poll(self,finished):
        try:
            while True:
                kind,data=self.queue.get_nowait()
                if kind=='progress':self.progress.configure(value=data[0]);self.message.set(data[1]);continue
                self.job=None;self.refresh();self.progress.configure(value=100 if kind=='done' else 0)
                if kind=='error':self.fail(data)
                else:
                    try:finished(data)
                    except Exception as exc:self.fail(str(exc))
                return
        except Empty:self.root.after(60,self.poll,finished)

    def fail(self,message):
        self.message.set(message);messagebox.showwarning('AirMem',message,parent=self.root)

    def close(self):
        if self.job is not None:self.message.set('Wait for the current operation to finish before closing.')
        else:self.root.destroy()


def start():
    root=tk.Tk();Window(root);root.mainloop();return 0
