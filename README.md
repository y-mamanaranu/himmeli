# Himmeli
Plot Himmeli in 3D &amp; 2D to create  mobiles. 

## Usage

```
him = Cone(60, 80, 3, folder=Path("Himmeli"))

fig, ax = plt.subplots(subplot_kw=dict(projection="3d"))
him.plot(ax)
if save:
    plt.savefig(him.file_3d)
plt.show()

fig, ax = him.plot_paper()
if save:
    plt.savefig(him.file_2d)
plt.show()
```

![Cone-3D](demo/Cone-3D.png)
![Cone-2D](demo/Cone-2D.png)

```
him = Bicone(45, 30, 35, 5, folder=Path("Himmeli"))

fig, ax = plt.subplots(subplot_kw=dict(projection="3d"))
him.plot(ax)
if save:
    plt.savefig(him.file_3d)
plt.show()

fig, ax = him.plot_paper()
if save:
    plt.savefig(him.file_2d)
plt.show()
```

![Bicone-3D](demo/Bicone-3D.png)
![Bicone-2D](demo/Bicone-2D.png)

```
him = Bicone_Lack(60, 30, 60, 20, 4, folder=Path("Himmeli"))

fig, ax = plt.subplots(subplot_kw=dict(projection="3d"))
him.plot(ax)
if save:
    plt.savefig(him.file_3d)
plt.show()

fig, ax = him.plot_paper()
if save:
    plt.savefig(him.file_2d)
plt.show()
```

![Bicone_Lack-3D](demo/Bicone_Lack-3D.png)
![Bicone_Lack-2D](demo/Bicone_Lack-2D.png)