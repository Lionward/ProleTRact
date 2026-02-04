import React, { useState, useEffect, useCallback } from 'react';
import axios from 'axios';
import './FileBrowser.css';

const API_BASE = process.env.REACT_APP_API_URL || 'http://localhost:8502';

export type FileBrowserMode = 'file' | 'folder';

interface BrowseEntry {
  name: string;
  path: string;
}

interface BrowseResponse {
  path: string;
  parent: string | null;
  directories: BrowseEntry[];
  files: BrowseEntry[];
}

interface FileBrowserProps {
  isOpen: boolean;
  onClose: () => void;
  onSelect: (path: string) => void;
  mode: FileBrowserMode;
  initialPath?: string;
  title?: string;
}

const FileBrowser: React.FC<FileBrowserProps> = ({
  isOpen,
  onClose,
  onSelect,
  mode,
  initialPath = '',
  title,
}) => {
  const [currentPath, setCurrentPath] = useState('');
  const [directories, setDirectories] = useState<BrowseEntry[]>([]);
  const [files, setFiles] = useState<BrowseEntry[]>([]);
  const [parent, setParent] = useState<string | null>(null);
  const [loading, setLoading] = useState(false);
  const [error, setError] = useState<string | null>(null);
  const [searchQuery, setSearchQuery] = useState('');

  const fetchContents = useCallback(async (path: string) => {
    setLoading(true);
    setError(null);
    try {
      const params = new URLSearchParams();
      if (path) params.set('path', path);
      params.set('mode', mode);
      const response = await axios.get<BrowseResponse>(
        `${API_BASE}/api/files/browse?${params.toString()}`
      );
      setCurrentPath(response.data.path);
      setDirectories(response.data.directories);
      setFiles(response.data.files);
      setParent(response.data.parent);
    } catch (err: unknown) {
      const msg = axios.isAxiosError(err)
        ? err.response?.data?.detail || err.message
        : String(err);
      setError(msg);
      setDirectories([]);
      setFiles([]);
    } finally {
      setLoading(false);
    }
  }, [mode]);

  useEffect(() => {
    if (isOpen) {
      setSearchQuery('');
      fetchContents(initialPath);
    }
  }, [isOpen, initialPath, fetchContents]);

  const q = searchQuery.trim().toLowerCase();
  const filteredDirs = q
    ? directories.filter((d) => d.name.toLowerCase().includes(q))
    : directories;
  const filteredFiles = q
    ? files.filter((f) => f.name.toLowerCase().includes(q))
    : files;

  const handleNavigate = (path: string) => {
    setSearchQuery('');
    fetchContents(path);
  };

  const handleSelectFile = (path: string) => {
    onSelect(path);
    onClose();
  };

  const handleSelectFolder = () => {
    if (currentPath) {
      onSelect(currentPath);
      onClose();
    }
  };

  const handleBackdropClick = (e: React.MouseEvent) => {
    if (e.target === e.currentTarget) onClose();
  };

  const handleKeyDown = (e: React.KeyboardEvent) => {
    if (e.key === 'Escape') onClose();
  };

  if (!isOpen) return null;

  const displayTitle = title ?? (mode === 'file' ? 'Select VCF File' : 'Select Folder');

  return (
    <div
      className="file-browser-overlay"
      onClick={handleBackdropClick}
      onKeyDown={handleKeyDown}
      role="dialog"
      aria-modal="true"
      aria-label={displayTitle}
    >
      <div className="file-browser-modal">
        <div className="file-browser-header">
          <h2 className="file-browser-title">{displayTitle}</h2>
          <button
            type="button"
            className="file-browser-close"
            onClick={onClose}
            aria-label="Close"
          >
            ×
          </button>
        </div>

        {/* Search */}
        <div className="file-browser-search">
          <input
            type="text"
            value={searchQuery}
            onChange={(e) => setSearchQuery(e.target.value)}
            placeholder="Search by name..."
            className="file-browser-search-input"
            autoComplete="off"
            autoFocus
          />
        </div>

        {/* Breadcrumb */}
        <div className="file-browser-breadcrumb">
          <button
            type="button"
            className="file-browser-breadcrumb-item"
            onClick={() => handleNavigate('')}
          >
            📁 Roots
          </button>
          {currentPath && (
            <>
              <span className="file-browser-breadcrumb-sep">/</span>
              <span className="file-browser-breadcrumb-current" title={currentPath}>
                {currentPath.split('/').filter(Boolean).slice(-2).join('/') || currentPath}
              </span>
            </>
          )}
        </div>

        {/* Content */}
        <div className="file-browser-content">
          {loading && (
            <div className="file-browser-loading">
              <span className="spinner" />
              <span>Loading...</span>
            </div>
          )}
          {error && (
            <div className="file-browser-error">
              {error}
            </div>
          )}
          {!loading && !error && (
            <div className="file-browser-list">
              {parent && parent !== currentPath && (
                <button
                  type="button"
                  className="file-browser-item file-browser-item-dir"
                  onClick={() => handleNavigate(parent)}
                >
                  <span className="file-browser-icon">📂</span>
                  <span>..</span>
                </button>
              )}
              {filteredDirs.map((d) => (
                <button
                  key={d.path}
                  type="button"
                  className="file-browser-item file-browser-item-dir"
                  onClick={() => handleNavigate(d.path)}
                >
                  <span className="file-browser-icon">📁</span>
                  <span>{d.name}</span>
                </button>
              ))}
              {filteredFiles.map((f) =>
                mode === 'file' ? (
                  <button
                    key={f.path}
                    type="button"
                    className="file-browser-item file-browser-item-file"
                    onClick={() => handleSelectFile(f.path)}
                  >
                    <span className="file-browser-icon">📄</span>
                    <span>{f.name}</span>
                  </button>
                ) : (
                  <div
                    key={f.path}
                    className="file-browser-item file-browser-item-file file-browser-item-readonly"
                  >
                    <span className="file-browser-icon">📄</span>
                    <span>{f.name}</span>
                  </div>
                )
              )}
              {!loading && !error && filteredDirs.length === 0 && filteredFiles.length === 0 && (
                <div className="file-browser-empty">
                  {q ? `No matches for "${searchQuery}"` : 'No items'}
                </div>
              )}
            </div>
          )}
        </div>

        {/* Footer */}
        <div className="file-browser-footer">
          {mode === 'folder' && (
            <button
              type="button"
              className="btn btn-primary btn-modern"
              onClick={handleSelectFolder}
              disabled={!currentPath}
            >
              Select this folder
            </button>
          )}
          <button type="button" className="btn btn-secondary btn-modern" onClick={onClose}>
            Cancel
          </button>
        </div>
      </div>
    </div>
  );
};

export default FileBrowser;
